#include "../third_party/htslib/htslib/bgzf.h"
#include "codec.hpp"
#include "io.hpp"
#include <algorithm>
#include <fstream>
#include <queue>
#include <stdexcept>
namespace nustar::star2711b {
namespace {
using bam::Bytes;
struct Record {
    std::uint64_t read = 0, order = 0;
    Bytes bytes;
    auto key() const {
        return std::tuple{bam::read32(bytes.data() + 4), bam::read32(bytes.data() + 8), read,
                          order};
    }
};
void writeRecord(std::ostream& out, const Record& r) {
    out.write(reinterpret_cast<const char*>(&r.read), sizeof(r.read));
    out.write(reinterpret_cast<const char*>(&r.order), sizeof(r.order));
    out.write(reinterpret_cast<const char*>(r.bytes.data()), r.bytes.size());
    if (!out)
        throw std::runtime_error("cannot write BAM sort temporary file");
}
bool readRecord(std::istream& in, Record& r) {
    if (!in.read(reinterpret_cast<char*>(&r.read), sizeof(r.read))) {
        if (in.eof() && in.gcount() == 0)
            return false;
        throw std::runtime_error("truncated BAM sort temporary file");
    }
    in.read(reinterpret_cast<char*>(&r.order), sizeof(r.order));
    r.bytes.resize(4);
    in.read(reinterpret_cast<char*>(r.bytes.data()), 4);
    const auto n = bam::read32(r.bytes.data());
    r.bytes.resize(n + 4);
    in.read(reinterpret_cast<char*>(r.bytes.data() + 4), n);
    if (!in)
        throw std::runtime_error("truncated BAM sort temporary record");
    return true;
}
struct BgzfCloser {
    void operator()(BGZF* p) const {
        if (p)
            bgzf_close(p);
    }
};
using Bgzf = std::unique_ptr<BGZF, BgzfCloser>;
void write(BGZF* file, std::span<const std::byte> bytes) {
    if (bgzf_write(file, bytes.data(), bytes.size()) != std::int64_t(bytes.size()))
        throw std::runtime_error("cannot write BAM output");
}
void close(Bgzf& file) {
    if (file && bgzf_close(file.release()) != 0)
        throw std::runtime_error("cannot finish BAM output");
}
struct Pending final : StagedOutput {
    OutputBatch output;
    explicit Pending(OutputBatch&& b) : StagedOutput(b.batch), output(std::move(b)) {}
};
class Sink final : public IOutputSink {
    OutputFiles files_;
    BamConfig config_;
    const ReferenceMetadata& ref_;
    bam::Encoder encoder_;
    BamCorrection correction_;
    std::unique_ptr<IOutputSink> auxiliary_;
    std::string prefix_;
    std::filesystem::path temporary_;
    Bgzf unsorted_;
    Bytes header_;
    std::vector<Record> records_;
    std::uint64_t memory_ = 0, serial_ = 0;
    std::vector<std::filesystem::path> runs_;
    Bgzf open(const std::string& path) {
        Bgzf file(bgzf_open(
            path.c_str(),
            (config_.compression < 0 ? std::string("w") : "w" + std::to_string(config_.compression))
                .c_str()));
        if (!file)
            throw std::runtime_error("cannot open BAM output: " + path);
        if (config_.threads > 1 && bgzf_mt(file.get(), config_.threads, 64) != 0)
            throw std::runtime_error("cannot start BAM compression workers");
        return file;
    }
    void flushRun() {
        if (records_.empty())
            return;
        std::sort(records_.begin(), records_.end(),
                  [](const auto& a, const auto& b) { return a.key() < b.key(); });
        const auto path = temporary_ / std::to_string(serial_++);
        std::ofstream file(path, std::ios::binary);
        for (const auto& r : records_)
            writeRecord(file, r);
        file.close();
        if (!file)
            throw std::runtime_error("cannot close BAM sort temporary file");
        runs_.push_back(path);
        records_.clear();
        memory_ = 0;
    }
    template <class Emit> void merge(std::span<const std::filesystem::path> paths, Emit emit) {
        std::vector<std::ifstream> files;
        files.reserve(paths.size());
        std::vector<Record> heads(paths.size());
        auto later = [&](std::size_t a, std::size_t b) { return heads[a].key() > heads[b].key(); };
        std::priority_queue<std::size_t, std::vector<std::size_t>, decltype(later)> queue(later);
        for (std::size_t i = 0; i < paths.size(); ++i) {
            files.emplace_back(paths[i], std::ios::binary);
            if (!files.back())
                throw std::runtime_error("cannot open BAM sort temporary file");
            if (readRecord(files.back(), heads[i]))
                queue.push(i);
        }
        while (!queue.empty()) {
            const auto i = queue.top();
            queue.pop();
            emit(heads[i]);
            if (readRecord(files[i], heads[i]))
                queue.push(i);
        }
    }

  public:
    Sink(const OutputFiles& files, const BamConfig& config, const ReferenceMetadata& ref,
         BamCorrection correction)
        : files_(files), config_(config), ref_(ref), encoder_(ref),
          correction_(std::move(correction)),
          prefix_(files.outFileNamePrefix.value_or(files.sam.string() + ".")),
          temporary_(prefix_ + "_NuSTARbamtmp") {
        auto auxfiles = files_;
        auxfiles.sam_enabled = false;
        auto result = bufferedSamSinkBackend().open(auxfiles);
        if (auto* error = std::get_if<Failure>(&result))
            throw std::runtime_error(error->message);
        auxiliary_ = std::get<std::unique_ptr<IOutputSink>>(std::move(result));
        if (config_.coordinate) {
            if (!std::filesystem::create_directory(temporary_))
                throw std::runtime_error("BAM temporary directory already exists: " +
                                         temporary_.string());
        }
    }
    ~Sink() override {
        std::error_code ignored;
        if (config_.coordinate)
            std::filesystem::remove_all(temporary_, ignored);
    }
    Result<EncodedFragment> encode(const FragmentId& id,
                                   const ProjectedOutput& output) const override {
        try {
            EncodedFragment encoded{id, {}, {}};
            std::size_t estimate = 0;
            for (const auto& record : output.records)
                estimate += 36 + record.qname.size() + 1 + 4 * record.cigar.size() +
                            (3 * record.sequence.size() + 1) / 2 + 128;
            encoded.sam.reserve(estimate);
            for (const auto& record : output.records)
                encoder_.append(encoded.sam, record);
            return encoded;
        } catch (const std::exception& e) {
            return Failure{ErrorKind::Io, Stage::Project, id, 1, e.what()};
        }
    }
    Status begin(std::span<const std::byte> sam_header) override {
        try {
            const auto status = auxiliary_->begin(sam_header);
            if (std::holds_alternative<Failure>(status))
                return status;
            std::string text(reinterpret_cast<const char*>(sam_header.data()), sam_header.size());
            // STAR emits separate HD sort orders for the two output files.
            auto makeHeader = [&](bool sorted) {
                std::string h = text;
                if (sorted) {
                    const auto end = h.find('\n');
                    h.replace(0, end, "@HD\tVN:1.4\tSO:coordinate");
                }
                Bytes bytes{std::byte('B'), std::byte('A'), std::byte('M'), std::byte(1)};
                bam::append32(bytes, h.size());
                const auto* p = reinterpret_cast<const std::byte*>(h.data());
                bytes.insert(bytes.end(), p, p + h.size());
                bam::append32(bytes, ref_.chromosomes.size());
                for (const auto& chr : ref_.chromosomes) {
                    bam::append32(bytes, chr.name.size() + 1);
                    p = reinterpret_cast<const std::byte*>(chr.name.c_str());
                    bytes.insert(bytes.end(), p, p + chr.name.size() + 1);
                    bam::append32(bytes, chr.length);
                }
                return bytes;
            };
            header_ = makeHeader(config_.coordinate);
            if (config_.unsorted) {
                unsorted_ = open(prefix_ + "Aligned.out.bam");
                write(unsorted_.get(), makeHeader(false));
            }
            return std::monostate{};
        } catch (const std::exception& e) {
            return Failure{ErrorKind::Io, Stage::Commit, {}, 1, e.what()};
        }
    }
    bool supportsStaging() const noexcept override { return true; }
    Result<std::unique_ptr<StagedOutput>> stage(OutputBatch&& output) override {
        return std::unique_ptr<StagedOutput>(std::make_unique<Pending>(std::move(output)));
    }
    Status submitStaged(std::unique_ptr<StagedOutput> p) override {
        return submit(std::move(static_cast<Pending&>(*p).output));
    }
    Status submit(OutputBatch&& output) override {
        try {
            for (const auto& fragment : output.fragments) {
                if (unsorted_)
                    write(unsorted_.get(), fragment.sam);
                if (!config_.coordinate)
                    continue;
                std::size_t offset = 0;
                std::uint64_t order = 0;
                while (offset < fragment.sam.size()) {
                    const auto n = bam::read32(fragment.sam.data() + offset) + 4;
                    Record r{fragment.fragment.ordinal,
                             order++,
                             {fragment.sam.begin() + offset, fragment.sam.begin() + offset + n}};
                    memory_ += n + sizeof(Record) * 2;
                    records_.push_back(std::move(r));
                    offset += n;
                    // Keep sorting records from inflating the mapping peak. The
                    // user budget is an upper bound, not a reservation.
                    if (memory_ >= std::min<std::uint64_t>(config_.sort_memory,
                                                           (8ULL << 20) * config_.threads))
                        flushRun();
                }
            }
            return std::monostate{};
        } catch (const std::exception& e) {
            return Failure{ErrorKind::Io, Stage::Commit, {}, 1, e.what()};
        }
    }
    Status finish(std::span<const std::byte> sj, const WorkerStatistics& statistics) override {
        try {
            close(unsorted_);
            return auxiliary_->finish(sj, statistics);
        } catch (const std::exception& e) {
            return Failure{ErrorKind::Io, Stage::Finish, {}, 1, e.what()};
        }
    }
    Status finalize() override {
        try {
            if (!config_.coordinate)
                return std::monostate{};
            auto sorted = open(prefix_ + "Aligned.sortedByCoord.out.bam");
            write(sorted.get(), header_);
            auto emit = [&](Record& record) {
                if (correction_) {
                    auto [cb, umi] = correction_(record.read);
                    bam::appendStringTag(record.bytes, "CB", cb);
                    bam::appendStringTag(record.bytes, "UB", umi);
                    const auto n = record.bytes.size() - 4;
                    for (unsigned i = 0; i < 4; ++i)
                        record.bytes[i] = std::byte(n >> (8 * i));
                }
                write(sorted.get(), record.bytes);
            };
            if (runs_.empty()) {
                std::sort(records_.begin(), records_.end(),
                          [](const auto& a, const auto& b) { return a.key() < b.key(); });
                for (auto& r : records_)
                    emit(r);
            } else {
                flushRun();
                // Bounded fan-in also bounds open descriptors for arbitrarily large input.
                while (runs_.size() > 64) {
                    std::vector<std::filesystem::path> next;
                    for (std::size_t i = 0; i < runs_.size(); i += 64) {
                        auto paths = std::span(runs_).subspan(
                            i, std::min<std::size_t>(64, runs_.size() - i));
                        const auto path = temporary_ / std::to_string(serial_++);
                        std::ofstream out(path, std::ios::binary);
                        merge(paths, [&](Record& r) { writeRecord(out, r); });
                        out.close();
                        if (!out)
                            throw std::runtime_error("cannot finish BAM merge temporary file");
                        next.push_back(path);
                        for (const auto& p : paths)
                            std::filesystem::remove(p);
                    }
                    runs_ = std::move(next);
                }
                merge(runs_, emit);
            }
            close(sorted);
            records_.clear();
            return std::monostate{};
        } catch (const std::exception& e) {
            return Failure{ErrorKind::Io, Stage::Finish, {}, 1, e.what()};
        }
    }
    void cancel() noexcept override { auxiliary_->cancel(); }
};
} // namespace
Result<std::unique_ptr<IOutputSink>> openBamSink(const OutputFiles& files, const BamConfig& config,
                                                 const ReferenceMetadata& reference,
                                                 BamCorrection correction) {
    try {
        return std::unique_ptr<IOutputSink>(
            std::make_unique<Sink>(files, config, reference, std::move(correction)));
    } catch (const std::exception& e) {
        return Failure{ErrorKind::Io, Stage::Commit, {}, 1, e.what()};
    }
}
} // namespace nustar::star2711b
