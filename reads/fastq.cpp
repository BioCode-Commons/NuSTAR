#include "raw.hpp"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

namespace nustar::star2711b {
namespace {

Failure inputFailure(std::string message) {
    return Failure{ErrorKind::Input, Stage::Input, {}, 1, std::move(message)};
}

class LineReader {
public:
    virtual ~LineReader() = default;
    virtual bool line(std::string& result) = 0;
};

class MappedLines final : public LineReader {
public:
    explicit MappedLines(const std::filesystem::path& path) {
        const auto filename = path.string();
        const int fd = ::open(filename.c_str(), O_RDONLY);
        if (fd < 0) throw std::runtime_error("cannot open FASTQ: " + filename);
        struct stat st {};
        if (::fstat(fd, &st) != 0 || !S_ISREG(st.st_mode) || st.st_size < 0) {
            ::close(fd);
            throw std::runtime_error("mmap requires a regular FASTQ file: " + filename);
        }
        size_ = static_cast<std::size_t>(st.st_size);
        const auto page_size = ::sysconf(_SC_PAGESIZE);
        if (page_size > 0) page_size_ = static_cast<std::size_t>(page_size);
        if (size_) {
            data_ = static_cast<const char*>(::mmap(nullptr, size_, PROT_READ,
                                                     MAP_PRIVATE, fd, 0));
            if (data_ == MAP_FAILED) {
                data_ = nullptr;
                ::close(fd);
                throw std::runtime_error("cannot mmap FASTQ: " + filename);
            }
        }
        ::close(fd);
    }
    ~MappedLines() override { if (data_) ::munmap(const_cast<char*>(data_), size_); }
    bool line(std::string& result) override {
        if (offset_ == size_) return false;
        const auto* begin = data_ + offset_;
        const auto* end = static_cast<const char*>(
            std::memchr(begin, '\n', size_ - offset_));
        const auto length = end ? static_cast<std::size_t>(end-begin) : size_-offset_;
        result.assign(begin, length && begin[length-1] == '\r' ? length-1 : length);
        offset_ += length + (end != nullptr);
        // line() returns an owning string; no consumer retains a pointer into
        // these completed pages. Drop old PTEs in large batches so mapping a
        // 200 GB FASTQ does not retain all consumed pages in the process RSS.
        // This is advisory and does not evict the shared filesystem cache.
        if (page_size_ && offset_-released_ >= (16ULL<<20)) {
            const auto consumed = offset_/page_size_*page_size_;
            if (::madvise(const_cast<char*>(data_)+released_, consumed-released_,
                          MADV_DONTNEED) == 0)
                released_ = consumed;
            else page_size_ = 0; // A failed optimization must not change input semantics.
        }
        return true;
    }
private:
    const char* data_ = nullptr;
    std::size_t size_ = 0, offset_ = 0;
    std::size_t page_size_ = 0, released_ = 0;
};

class BufferedLines final : public LineReader {
public:
    explicit BufferedLines(const std::filesystem::path& path) : input_(path) {
        if (!input_) throw std::runtime_error("cannot open FASTQ: " + path.string());
    }
    bool line(std::string& result) override {
        if (!std::getline(input_, result)) {
            if (input_.bad()) throw std::runtime_error("failed reading FASTQ");
            return false;
        }
        if (!result.empty() && result.back() == '\r') result.pop_back();
        return true;
    }
private:
    std::ifstream input_;
};

using OpenLines = std::unique_ptr<LineReader>(*)(const std::filesystem::path&);

template<class T> std::unique_ptr<LineReader> createLines(const std::filesystem::path& path) {
    return std::make_unique<T>(path);
}

class FastqReadSource final : public IReadSource {
public:
    FastqReadSource(InputSpec spec, OpenLines opener)
        : spec_(std::move(spec)), opener_(opener) {
        if (spec_.groups.empty()) throw std::runtime_error("no FASTQ input groups");
        if (!spec_.readFilesCommand.empty())
            throw std::runtime_error("readFilesCommand is not supported by this FASTQ backend");
        openGroup();
    }
    Result<std::optional<ReadBatch>> next(std::uint64_t target_bytes) override {
        try {
            if (cancelled_) return inputFailure("FASTQ input cancelled");
            if (group_ == spec_.groups.size()) return std::optional<ReadBatch>{};
            if (spec_.readMapNumber && ordinal_ >= *spec_.readMapNumber)
                return std::optional<ReadBatch>{};
            ReadBatch batch{BatchId{next_batch_}, {}};
            std::uint64_t bytes = 0;
            const auto target = std::max<std::uint64_t>(1, target_bytes);
            while (bytes < target) {
                if (group_ == spec_.groups.size()) break;
                if (spec_.readMapNumber && ordinal_ >= *spec_.readMapNumber) break;
                RawMate first;
                if (!readRecord(*first_, first)) {
                    if (spec_.layout == Layout::PE) {
                        RawMate excess;
                        if (readRecord(*second_, excess))
                            throw std::runtime_error("mate 2 has extra FASTQ records");
                    }
                    ++group_;
                    if (group_ == spec_.groups.size()) break;
                    openGroup();
                    continue;
                }
                RawFragment fragment{FragmentId{static_cast<std::uint32_t>(group_), ordinal_++},
                    spec_.layout, FastxFormat::Fastq, {}, 'N'};
                bytes += first.header.size()+first.sequence.size()+first.quality->size()+2;
                fragment.mates[0] = std::move(first);
                if (spec_.layout == Layout::PE) {
                    RawMate second;
                    if (!readRecord(*second_, second))
                        throw std::runtime_error("mate 2 ended before mate 1");
                    bytes += second.header.size()+second.sequence.size()+second.quality->size()+2;
                    fragment.mates[1] = std::move(second);
                }
                batch.fragments.push_back(std::move(fragment));
            }
            if (batch.fragments.empty()) return std::optional<ReadBatch>{};
            ++next_batch_;
            return std::optional<ReadBatch>{std::move(batch)};
        } catch (const std::exception& error) {
            return inputFailure(error.what());
        }
    }
    void cancel() noexcept override { cancelled_ = true; }
private:
    static bool readRecord(LineReader& lines, RawMate& mate) {
        if (!lines.line(mate.header)) return false;
        std::string plus;
        std::string quality;
        // STAR consumes the third line without validating its contents.
        if (!lines.line(mate.sequence) || !lines.line(plus) || !lines.line(quality) ||
            mate.header.empty() || mate.header.front() != '@' ||
            mate.sequence.size() != quality.size())
            throw std::runtime_error("invalid FASTQ record");
        mate.quality = std::move(quality);
        return true;
    }
    void openGroup() {
        first_ = opener_(spec_.groups[group_].mate_paths[0]);
        second_ = spec_.layout == Layout::PE ?
            opener_(spec_.groups[group_].mate_paths[1]) : nullptr;
    }
    InputSpec spec_;
    OpenLines opener_;
    std::unique_ptr<LineReader> first_, second_;
    std::size_t group_ = 0;
    std::uint64_t ordinal_ = 0, next_batch_ = 0;
    bool cancelled_ = false;
};

class FastqBackend final : public IReadSourceBackend {
public:
    explicit FastqBackend(OpenLines opener) : opener_(opener) {}
    Result<std::unique_ptr<IReadSource>> open(const InputSpec& spec) const override {
        try { return std::make_unique<FastqReadSource>(spec, opener_); }
        catch (const std::exception& error) { return inputFailure(error.what()); }
    }
private:
    OpenLines opener_;
};

} // namespace

const IReadSourceBackend& memoryMappedFastqBackend() {
    static const FastqBackend backend(&createLines<MappedLines>);
    return backend;
}
const IReadSourceBackend& bufferedFastqBackend() {
    static const FastqBackend backend(&createLines<BufferedLines>);
    return backend;
}
Result<std::unique_ptr<IReadSource>> openReadSource(const InputSpec& spec) {
    return memoryMappedFastqBackend().open(spec);
}

} // namespace nustar::star2711b
