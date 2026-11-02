#include "sink.hpp"

#include <fstream>
#include <array>
#include <atomic>
#include <mutex>

namespace nustar::star2711b {
namespace {

Failure ioFailure(std::string message) {
    return Failure{ErrorKind::Io, Stage::Commit, {}, 1, std::move(message)};
}

struct ChunkDirectory {
    std::filesystem::path path;
    explicit ChunkDirectory(const std::filesystem::path& sam) {
        static std::atomic<std::uint64_t> next{0};
        for (;;) {
            auto candidate=std::filesystem::path(sam.string()+".chunks."+std::to_string(next.fetch_add(1)));
            std::error_code ec;
            if (std::filesystem::create_directory(candidate,ec)) { path=std::move(candidate); break; }
            if (ec) throw std::runtime_error("cannot create SAM chunk directory: "+ec.message());
        }
    }
    ~ChunkDirectory() { std::error_code ec; std::filesystem::remove_all(path,ec); }
};
struct FileChunk final : StagedOutput {
    FileChunk(BatchId id,std::shared_ptr<ChunkDirectory> owner)
        : StagedOutput(id), directory(std::move(owner)),
          path(directory->path/std::to_string(id.value)) {}
    ~FileChunk() { std::error_code ec; std::filesystem::remove(path,ec); }
    std::shared_ptr<ChunkDirectory> directory;
    std::filesystem::path path;
};

class FileSamSink final : public IOutputSink {
public:
    explicit FileSamSink(OutputFiles files)
        : files_(std::move(files)) {
        if(files_.sam_enabled) sam_.open(files_.sam,std::ios::binary);
        if (files_.sam_enabled && !sam_) throw std::runtime_error("cannot open SAM output");
    }
    Status begin(std::span<const std::byte> header) override {
        if (begun_) return ioFailure("SAM sink already started");
        begun_ = true;
        return write(header);
    }
    Status submit(OutputBatch&& batch) override {
        if (!begun_ || batch.batch.value != next_batch_)
            return ioFailure("SAM batches must be submitted in input order");
        for (const auto& fragment : batch.fragments) {
            const auto status = write(fragment.sam);
            if (std::holds_alternative<Failure>(status)) return status;
        }
        ++next_batch_;
        return std::monostate{};
    }
    bool supportsStaging() const noexcept override { return files_.sam_enabled; }
    Result<std::unique_ptr<StagedOutput>> stage(OutputBatch&& batch) override {
        if (!begun_ || finished_) return ioFailure("SAM sink is not active");
        try {
            std::call_once(chunk_once_,[this] { chunks_=std::make_shared<ChunkDirectory>(files_.sam); });
            auto chunk=std::make_unique<FileChunk>(batch.batch,chunks_);
            std::array<char,1<<16> buffer;
            std::ofstream file;
            file.rdbuf()->pubsetbuf(buffer.data(),buffer.size());
            file.open(chunk->path,std::ios::binary);
            if (!file) return ioFailure("cannot open SAM chunk");
            for (const auto& fragment : batch.fragments)
                file.write(reinterpret_cast<const char*>(fragment.sam.data()),fragment.sam.size());
            file.close();
            if (!file) return ioFailure("failed to write SAM chunk");
            // Release the small per-read allocations on the producing worker.
            batch.fragments.clear();
            return std::unique_ptr<StagedOutput>(std::move(chunk));
        } catch (const std::exception& e) { return ioFailure(e.what()); }
    }
    Status submitStaged(std::unique_ptr<StagedOutput> staged) override {
        const auto* chunk=dynamic_cast<const FileChunk*>(staged.get());
        if (!begun_ || finished_ || !chunk || chunk->directory!=chunks_ ||
            chunk->batch.value!=next_batch_)
            return ioFailure("invalid or unordered SAM chunk");
        std::ifstream input(chunk->path,std::ios::binary);
        if (!input) return ioFailure("cannot read SAM chunk");
        std::array<char,1<<18> buffer;
        while (input) {
            input.read(buffer.data(),buffer.size());
            if (input.gcount()) sam_.write(buffer.data(),input.gcount());
        }
        if (!input.eof() || !sam_) return ioFailure("failed to concatenate SAM chunk");
        ++next_batch_;
        return std::monostate{};
    }
    Status finish(std::span<const std::byte> sj_table,
                  const WorkerStatistics& stats) override {
        if (!begun_ || finished_) return ioFailure("SAM sink is not active");
        if(files_.sam_enabled) {
            sam_.flush();
            if (!sam_) return ioFailure("failed to write SAM");
        }
        const auto finished = std::chrono::system_clock::now();
        std::ofstream sj(files_.outFileNamePrefix.value_or(files_.sam.string()+".")+"SJ.out.tab", std::ios::binary);
        if (!sj) return ioFailure("cannot open SJ output");
        sj.write(reinterpret_cast<const char*>(sj_table.data()),
                 static_cast<std::streamsize>(sj_table.size()));
        sj.close();
        if (!sj) return ioFailure("failed to write SJ table");
        std::ofstream log(files_.outFileNamePrefix.value_or(files_.sam.string()+".")+"Log.final.out");
        if (!log) return ioFailure("cannot open statistics output");
        log << formatStarFinalStatistics(stats.totals, files_.job_started,
                                         files_.mapping_started, finished);
        log.close();
        if (!log) return ioFailure("failed to write statistics");
        finished_ = true;
        return std::monostate{};
    }
    void cancel() noexcept override { sam_.close(); }
private:
    Status write(std::span<const std::byte> bytes) {
        if(!files_.sam_enabled) return std::monostate{};
        sam_.write(reinterpret_cast<const char*>(bytes.data()),
                   static_cast<std::streamsize>(bytes.size()));
        if (!sam_) return ioFailure("failed to write SAM");
        return std::monostate{};
    }
    std::once_flag chunk_once_;
    std::shared_ptr<ChunkDirectory> chunks_;
    OutputFiles files_;
    std::ofstream sam_;
    std::uint64_t next_batch_ = 0;
    bool begun_ = false, finished_ = false;
};

class FileBackend final : public IOutputSinkBackend {
public:
    Result<std::unique_ptr<IOutputSink>> open(const OutputFiles& files) const override {
        try { return std::make_unique<FileSamSink>(files); }
        catch (const std::exception& error) { return ioFailure(error.what()); }
    }
};

} // namespace

const IOutputSinkBackend& bufferedSamSinkBackend() {
    static const FileBackend backend;
    return backend;
}

} // namespace nustar::star2711b
