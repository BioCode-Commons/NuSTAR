#include "raw.hpp"

#include <algorithm>
#include <atomic>
#include <bit>
#include <cerrno>
#include <cstring>
#include <stdexcept>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

namespace nustar::star2711b {
namespace {

constexpr std::uint64_t unknown = UINT64_MAX;
constexpr std::size_t block_bytes = 4ULL << 20;
constexpr unsigned window_blocks = 32; // Per mate, across all input groups.
constexpr unsigned read_concurrency = 4;
constexpr unsigned max_workers = 256;

Failure inputFailure(std::string message) {
    return {ErrorKind::Input, Stage::Input, {}, 1, std::move(message)};
}

struct BlockData {
    std::vector<char> bytes;
    // Offsets immediately after newlines; EOF supplies the last line ending
    // when a nonempty final line has no trailing newline.
    std::vector<std::uint32_t> ends;
};

struct Block {
    // Metadata stays alive until the source is destroyed. Only active blocks
    // own bytes/line offsets, so helpers can safely finish an old prefix update
    // even after a consumer has released that block's payload.
    std::atomic<bool> ready{false};
    std::atomic<std::uint64_t> after{unknown};
    std::uint32_t lines = 0, last_end = 0;
    bool ends_line = false;
    unsigned slot = 0;
    BlockData* data = nullptr;
};

class File {
public:
    explicit File(const std::filesystem::path& path) {
        fd = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
        if (fd < 0) throw std::runtime_error("cannot open FASTQ: " + path.string());
        struct stat st{};
        if (::fstat(fd, &st) || !S_ISREG(st.st_mode) || st.st_size < 0) {
            ::close(fd);
            throw std::runtime_error("parallel input requires a regular FASTQ file: " + path.string());
        }
        size = static_cast<std::uint64_t>(st.st_size);
        count = size / block_bytes + (size % block_bytes != 0);
        try { blocks = std::make_unique<Block[]>(count); }
        catch (...) { ::close(fd); throw; }
    }
    ~File() { ::close(fd); }
    File(const File&) = delete;
    File& operator=(const File&) = delete;

    void read(std::uint64_t id, BlockData& storage, unsigned slot) {
        auto* data = &storage;
        data->ends.clear();
        const auto begin = id * block_bytes;
        data->bytes.resize(std::min<std::uint64_t>(block_bytes, size - begin));
        std::size_t done = 0;
        while (done < data->bytes.size()) {
            const auto n = ::pread(fd, data->bytes.data() + done,
                                   data->bytes.size() - done, begin + done);
            if (n < 0 && errno == EINTR) continue;
            if (n <= 0) throw std::runtime_error("failed reading FASTQ range: " +
                std::string(n < 0 ? std::strerror(errno) : "unexpected EOF"));
            done += static_cast<std::size_t>(n);
        }
        const char* first = data->bytes.data();
        const char* cursor = first;
        const char* end = first + data->bytes.size();
        while (const auto* newline = static_cast<const char*>(std::memchr(cursor, '\n', end - cursor))) {
            cursor = newline + 1;
            data->ends.push_back(static_cast<std::uint32_t>(cursor - first));
        }
        if (id + 1 == count && end[-1] != '\n')
            data->ends.push_back(static_cast<std::uint32_t>(end - first));
        auto& block = blocks[id];
        block.lines = static_cast<std::uint32_t>(data->ends.size());
        block.last_end = data->ends.empty() ? 0 : data->ends.back();
        block.ends_line = end[-1] == '\n' || id + 1 == count;
        block.slot = slot;
        block.data = data;
        block.ready.store(true, std::memory_order_release);
    }

    void advance() {
        auto id = prefix.load(std::memory_order_acquire);
        while (id < count && blocks[id].ready.load(std::memory_order_acquire)) {
            const auto before = id ? blocks[id - 1].after.load(std::memory_order_acquire) : 0;
            blocks[id].after.store(before + blocks[id].lines, std::memory_order_release);
            prefix.compare_exchange_strong(id, id + 1, std::memory_order_acq_rel);
            id = prefix.load(std::memory_order_acquire);
        }
    }
    std::uint64_t availableLines() const {
        const auto n = prefix.load(std::memory_order_acquire);
        return n ? blocks[n - 1].after.load(std::memory_order_acquire) : 0;
    }
    std::pair<std::uint64_t, bool> progress() const {
        const auto n = prefix.load(std::memory_order_acquire);
        return {n ? blocks[n - 1].after.load(std::memory_order_acquire) : 0, n == count};
    }

    int fd = -1;
    std::uint64_t size = 0, count = 0;
    std::unique_ptr<Block[]> blocks;
    std::atomic<std::uint64_t> next{0}, prefix{0}, retired{0};
};

// Searches the sparse block prefix once per batch. Every subsequent line is
// obtained directly from the already-computed local newline offsets.
class LineCursor {
public:
    LineCursor(File& file, std::uint64_t preceding_lines) : file_(file) {
        if (!preceding_lines) return;
        std::uint64_t lo = 0, hi = file.prefix.load(std::memory_order_acquire);
        while (lo < hi) {
            const auto mid = lo + (hi - lo) / 2;
            if (file.blocks[mid].after.load(std::memory_order_acquire) < preceding_lines) lo = mid + 1;
            else hi = mid;
        }
        block_ = lo;
        const auto before = lo ? file.blocks[lo - 1].after.load(std::memory_order_acquire) : 0;
        line_ = preceding_lines - before;
        const auto& block = file.blocks[lo];
        // The previous block may already be retired when the batch starts
        // exactly at its end. Its last line offset is persistent metadata.
        position_ = lo * block_bytes + (line_ == block.lines ? block.last_end : block.data->ends[line_ - 1]);
    }
    void skip() { position_ = nextEnd(); }
    void line(std::string& out) {
        const auto end = nextEnd();
        out.resize(end - position_);
        auto position = position_;
        std::size_t written = 0;
        while (position < end) {
            const auto id = position / block_bytes;
            const auto offset = position % block_bytes;
            const auto& bytes = file_.blocks[id].data->bytes;
            const auto n = std::min<std::uint64_t>(end - position, bytes.size() - offset);
            std::memcpy(out.data() + written, bytes.data() + offset, n);
            position += n;
            written += n;
        }
        if (!out.empty() && out.back() == '\n') out.pop_back();
        if (!out.empty() && out.back() == '\r') out.pop_back();
        position_ = end;
    }
private:
    std::uint64_t nextEnd() {
        while (line_ == file_.blocks[block_].lines) { ++block_; line_ = 0; }
        return block_ * block_bytes + file_.blocks[block_].data->ends[line_++];
    }
    File& file_;
    std::uint64_t block_ = 0, line_ = 0, position_ = 0;
};

struct Group {
    explicit Group(InputGroup value) : input(std::move(value)) {}
    InputGroup input;
    std::atomic<unsigned> initialized{0};
    std::array<std::unique_ptr<File>, 2> files;
    std::atomic<std::uint64_t> batch_base{unknown}, ordinal_base{unknown};
    std::atomic<std::uint64_t> end_batch{unknown};
};

class ParallelFastqSource final : public IReadSource {
public:
    ParallelFastqSource(const InputSpec& spec, std::uint64_t target, std::uint64_t limit)
        : layout_(spec.layout), mates_(spec.layout == Layout::PE ? 2 : 1),
          batch_records_(std::clamp<std::uint64_t>(target / (mates_ * 384), 1, 4096)),
          limit_(spec.readMapNumber ? std::min(limit, *spec.readMapNumber) : limit) {
        static_assert(std::atomic<std::uint64_t>::is_always_lock_free);
        for (auto& lease : leases_) lease.store(unknown, std::memory_order_relaxed);
        for (const auto& input : spec.groups) groups_.push_back(std::make_unique<Group>(input));
        if (!groups_.empty()) {
            groups_[0]->batch_base.store(0);
            groups_[0]->ordinal_base.store(0);
        }
    }
    Status prepareParallel(std::uint64_t bytes, std::uint32_t,
                           std::uint64_t limit = UINT64_MAX) override {
        batch_records_ = std::clamp<std::uint64_t>(bytes / (mates_ * 384), 1, 4096);
        limit_ = std::min(limit_, limit);
        return std::monostate{};
    }
    bool concurrentReads() const noexcept override { return true; }

    Result<std::optional<ReadBatch>> next(std::uint64_t) override {
        Lease lease(*this);
        try {
            for (;;) {
                const auto changed = changed_.load(std::memory_order_acquire);
                if (cancelled_.load(std::memory_order_relaxed)) return inputFailure("FASTQ input cancelled");
                auto ticket = next_batch_.load();
                lease.value.store(ticket);
                reclaim();
                const auto group_id = current_group_.load(std::memory_order_acquire);
                if (group_id == groups_.size()) return std::optional<ReadBatch>{};
                auto& group = *groups_[group_id];
                const auto base = group.batch_base.load(std::memory_order_acquire);
                const auto ordinal = group.ordinal_base.load(std::memory_order_acquire);
                if (ticket < base) continue;
                const auto remaining = limit_ - ordinal;
                if (!remaining) {
                    group.end_batch.store(base, std::memory_order_release);
                    current_group_.store(groups_.size(), std::memory_order_release);
                    signal();
                    continue;
                }
                if (!initialize(group)) { changed_.wait(changed, std::memory_order_acquire); continue; }
                for (unsigned mate = 0; mate < mates_; ++mate) group.files[mate]->advance();
                const auto first_progress = group.files[0]->progress();
                auto available = first_progress.first / 4;
                bool complete = first_progress.second;
                if (mates_ == 2) {
                    const auto second_progress = group.files[1]->progress();
                    available = std::min(available, second_progress.first / 4);
                    complete &= second_progress.second;
                }
                available = std::min(available, remaining);
                const auto start = (ticket - base) * batch_records_;
                auto count = std::min(batch_records_, remaining - std::min(start, remaining));
                if (complete && available < remaining) {
                    for (unsigned mate = 0; mate < mates_; ++mate)
                        if (group.files[mate]->availableLines() % 4)
                            throw std::runtime_error("incomplete FASTQ record at EOF");
                    if (mates_ == 2 && group.files[0]->availableLines() != group.files[1]->availableLines())
                        throw std::runtime_error("paired FASTQ files have different record counts");
                    count = std::min(count, available - std::min(start, available));
                }
                if (count && start <= available && count <= available - start) {
                    if (!next_batch_.compare_exchange_strong(ticket, ticket + 1)) continue;
                    ReadBatch batch{BatchId{ticket}, {}};
                    batch.fragments.reserve(count);
                    LineCursor first(*group.files[0], start * 4);
                    std::optional<LineCursor> second;
                    if (mates_ == 2) second.emplace(*group.files[1], start * 4);
                    for (std::uint64_t i = 0; i < count; ++i) {
                        RawFragment fragment{FragmentId{static_cast<std::uint32_t>(group_id), ordinal + start + i},
                                             layout_, FastxFormat::Fastq, {}, 'N'};
                        auto& a = fragment.mates[0].emplace();
                        first.line(a.header); first.line(a.sequence); first.skip(); first.line(a.quality.emplace());
                        if (second) {
                            auto& b = fragment.mates[1].emplace();
                            second->skip(); b.header = a.header; // STAR takes both mate IDs from R1.
                            second->line(b.sequence); second->skip(); second->line(b.quality.emplace());
                        }
                        batch.fragments.push_back(std::move(fragment));
                    }
                    lease.release();
                    reclaim();
                    signal();
                    return std::optional<ReadBatch>{std::move(batch)};
                }
                if ((complete || available == remaining) && start >= available) {
                    const auto end = base + available / batch_records_ + (available % batch_records_ != 0);
                    group.end_batch.store(end, std::memory_order_release);
                    if (available == remaining || group_id + 1 == groups_.size()) {
                        auto expected = group_id;
                        current_group_.compare_exchange_strong(expected, groups_.size(), std::memory_order_acq_rel);
                    } else {
                        auto& following = *groups_[group_id + 1];
                        following.batch_base.store(end, std::memory_order_release);
                        following.ordinal_base.store(ordinal + available, std::memory_order_release);
                        auto expected = group_id;
                        current_group_.compare_exchange_strong(expected, group_id + 1, std::memory_order_acq_rel);
                    }
                    signal();
                    continue;
                }
                unsigned preferred = mates_ == 2 && group.files[1]->availableLines() < group.files[0]->availableLines();
                if (readBlock(group, preferred) || (mates_ == 2 && readBlock(group, 1 - preferred))) continue;
                // Fail resource exhaustion instead of waiting forever on a
                // batch whose lines cannot fit in the bounded input window.
                if (readers_.load() == 0 && !copiesInFlight()) {
                    for (unsigned mate = 0; mate < mates_; ++mate) {
                        auto& file = *group.files[mate];
                        if (buffered_[mate].load() == window_blocks &&
                            file.prefix.load() == file.next.load() &&
                            file.availableLines() / 4 < start + count)
                            throw std::runtime_error("FASTQ batch exceeds the parallel input buffer window");
                    }
                }
                changed_.wait(changed, std::memory_order_acquire);
            }
        } catch (const std::exception& error) {
            cancel();
            return inputFailure(error.what());
        }
    }
    void cancel() noexcept override {
        cancelled_.store(true, std::memory_order_relaxed);
        signal();
    }

private:
    // Announce a ticket before reserving it. Reclamation snapshots next_batch
    // first, then all leases (sequentially consistent), so payloads stay alive
    // until every caller copying their strings has finished.
    struct Lease {
        explicit Lease(ParallelFastqSource& source) : owner(source), value(acquire(source)) {}
        ~Lease() { release(); }
        void release() {
            if (active) { value.store(unknown); active = false; owner.signal(); }
        }
        static std::atomic<std::uint64_t>& acquire(ParallelFastqSource& source) {
            for (;;) {
                const auto changed = source.changed_.load();
                for (auto& slot : source.leases_) {
                    auto empty = unknown;
                    if (slot.compare_exchange_strong(empty, 0)) return slot;
                }
                source.changed_.wait(changed);
            }
        }
        ParallelFastqSource& owner;
        std::atomic<std::uint64_t>& value;
        bool active = true;
    };

    void signal() noexcept { changed_.fetch_add(1, std::memory_order_release); changed_.notify_all(); }
    bool initialize(Group& group) {
        auto state = group.initialized.load(std::memory_order_acquire);
        if (state == 2) return true;
        if (state == 0 && group.initialized.compare_exchange_strong(state, 1, std::memory_order_acq_rel)) {
            for (unsigned mate = 0; mate < mates_; ++mate)
                group.files[mate] = std::make_unique<File>(group.input.mate_paths[mate]);
            group.initialized.store(2, std::memory_order_release);
            signal();
            return true;
        }
        return false;
    }
    bool readBlock(Group& group, unsigned mate) {
        auto& file = *group.files[mate];
        if (file.next.load(std::memory_order_relaxed) >= file.count) return false;
        auto readers = readers_.load();
        while (readers < read_concurrency && !readers_.compare_exchange_weak(readers, readers + 1)) {}
        if (readers >= read_concurrency) return false;
        auto buffered = buffered_[mate].load();
        while (buffered < window_blocks && !buffered_[mate].compare_exchange_weak(buffered, buffered + 1)) {}
        if (buffered >= window_blocks) { readers_.fetch_sub(1); return false; }
        auto id = file.next.load();
        while (id < file.count && !file.next.compare_exchange_weak(id, id + 1)) {}
        if (id == file.count) {
            buffered_[mate].fetch_sub(1); readers_.fetch_sub(1); signal(); return false;
        }
        // Each reserved buffer credit guarantees a free pool slot. Recycling
        // keeps allocation/unmapping out of the steady-state read loop.
        auto free = free_slots_.load(std::memory_order_acquire);
        unsigned slot;
        do { slot = std::countr_zero(free); }
        while (!free_slots_.compare_exchange_weak(free, free & ~(1ULL << slot), std::memory_order_acq_rel));
        try { file.read(id, pool_[slot], slot); file.advance(); }
        catch (...) {
            free_slots_.fetch_or(1ULL << slot, std::memory_order_release);
            buffered_[mate].fetch_sub(1); readers_.fetch_sub(1); signal(); throw;
        }
        readers_.fetch_sub(1);
        signal();
        return true;
    }
    bool copiesInFlight() const {
        const auto next = next_batch_.load();
        for (const auto& lease : leases_) if (lease.load() < next) return true;
        return false;
    }
    void reclaim() {
        auto first_needed = next_batch_.load();
        for (const auto& lease : leases_) first_needed = std::min(first_needed, lease.load());
        const auto active = current_group_.load(std::memory_order_acquire);
        for (std::size_t index = 0; index < groups_.size() && index <= active; ++index) {
            auto& group = *groups_[index];
            if (group.initialized.load(std::memory_order_acquire) != 2) continue;
            const auto base = group.batch_base.load(std::memory_order_acquire);
            if (first_needed < base) continue;
            const bool finished = first_needed >= group.end_batch.load(std::memory_order_acquire);
            const auto lines = (first_needed - base) * batch_records_ * 4;
            for (unsigned mate = 0; mate < mates_; ++mate) {
                auto& file = *group.files[mate];
                auto id = file.retired.load();
                while (id < file.count && file.blocks[id].ready.load(std::memory_order_acquire)) {
                    auto& block = file.blocks[id];
                    const auto after = block.after.load(std::memory_order_acquire);
                    if (!finished && (after == unknown || after > lines || (after == lines && !block.ends_line))) break;
                    if (file.retired.compare_exchange_strong(id, id + 1)) {
                        block.data = nullptr;
                        free_slots_.fetch_or(1ULL << block.slot, std::memory_order_release);
                        buffered_[mate].fetch_sub(1);
                        ++id;
                        signal();
                    }
                }
            }
        }
    }

    Layout layout_;
    unsigned mates_;
    std::uint64_t batch_records_, limit_;
    std::vector<std::unique_ptr<Group>> groups_;
    alignas(64) std::atomic<std::uint64_t> next_batch_{0};
    std::atomic<std::size_t> current_group_{0};
    std::array<std::atomic<std::uint64_t>, max_workers> leases_;
    std::array<std::atomic<unsigned>, 2> buffered_{};
    std::array<BlockData, window_blocks * 2> pool_;
    std::atomic<std::uint64_t> free_slots_{UINT64_MAX};
    std::atomic<unsigned> readers_{0};
    std::atomic<std::uint64_t> changed_{0};
    std::atomic<bool> cancelled_{false};
};

} // namespace

// Retain this entry point for callers supplying an explicit batch size and limit.
// All byte reads happen on demand inside next().
Result<std::unique_ptr<IReadSource>> openIndexedFastqSource(
    const InputSpec& spec, std::uint64_t target_bytes, std::uint64_t limit) {
    try { return std::unique_ptr<IReadSource>(std::make_unique<ParallelFastqSource>(spec, target_bytes, limit)); }
    catch (const std::exception& error) { return inputFailure(error.what()); }
}

namespace {
class FastqBackend final : public IReadSourceBackend {
public:
    Result<std::unique_ptr<IReadSource>> open(const InputSpec& spec) const override {
        if (!spec.readFilesCommand.empty())
            return inputFailure("readFilesCommand is not supported by this FASTQ backend");
        return openIndexedFastqSource(spec, 1ULL << 20, UINT64_MAX);
    }
};
} // namespace

const IReadSourceBackend& fastqSourceBackend() {
    static const FastqBackend backend;
    return backend;
}

// Compatibility names; both use the single on-demand pread implementation.
const IReadSourceBackend& memoryMappedFastqBackend() { return fastqSourceBackend(); }
const IReadSourceBackend& bufferedFastqBackend() { return fastqSourceBackend(); }
Result<std::unique_ptr<IReadSource>> openReadSource(const InputSpec& spec) {
    return fastqSourceBackend().open(spec);
}

} // namespace nustar::star2711b
