#include "sink.hpp"

#include <fstream>
#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cstring>
#include <limits>
#include <fcntl.h>
#include <unistd.h>

namespace nustar::star2711b {
namespace {

Failure ioFailure(std::string message) {
    return Failure{ErrorKind::Io, Stage::Commit, {}, 1, std::move(message)};
}

class FileSamSink final : public IOutputSink {
public:
    explicit FileSamSink(OutputFiles files)
        : files_(std::move(files)) {
        static_assert(std::atomic<std::uint64_t>::is_always_lock_free);
        if (files_.sam_enabled) {
            fd_ = ::open(files_.sam.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC, 0666);
            if (fd_ < 0) throw std::runtime_error("cannot open SAM output");
        }
    }
    ~FileSamSink() override { if (fd_ >= 0) ::close(fd_); }
    Status begin(std::span<const std::byte> header) override {
        return write(header);
    }
    bool concurrentWrites() const noexcept override {
        return true;
    }
    Status submit(OutputBatch&& batch) override {
        if (!files_.sam_enabled) return std::monostate{};
        std::size_t size = 0;
        for (const auto& fragment : batch.fragments)
            size += fragment.sam.size();
        // Coalesce on the producing worker: one range reservation and usually
        // one syscall per batch, preserving every fragment's record adjacency.
        std::vector<std::byte> bytes;
        bytes.reserve(size);
        for (const auto& fragment : batch.fragments)
            bytes.insert(bytes.end(), fragment.sam.begin(), fragment.sam.end());
        return write(bytes);
    }
    Status finish(std::span<const std::byte> sj_table,
                  const WorkerStatistics& stats) override {
        if (fd_ >= 0) {
            const int fd = fd_;
            fd_ = -1;
            if (::close(fd) != 0) return ioFailure("failed to close SAM output");
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
        return std::monostate{};
    }
    void cancel() noexcept override {
        if (fd_ >= 0) { ::close(fd_); fd_ = -1; }
    }
private:
    Status write(std::span<const std::byte> bytes) {
        if(!files_.sam_enabled) return std::monostate{};
        const auto offset = output_offset_.fetch_add(bytes.size(), std::memory_order_relaxed);
        std::size_t written = 0;
        while (written < bytes.size()) {
            const auto length = std::min<std::size_t>(bytes.size() - written,
                std::numeric_limits<ssize_t>::max());
            const auto n = ::pwrite(fd_, bytes.data() + written, length,
                                    static_cast<off_t>(offset + written));
            if (n < 0 && errno == EINTR) continue;
            if (n <= 0) return ioFailure("failed to write SAM range: " +
                std::string(n < 0 ? std::strerror(errno) : "zero-length write"));
            written += static_cast<std::size_t>(n);
        }
        return std::monostate{};
    }
    OutputFiles files_;
    int fd_ = -1;
    alignas(64) std::atomic<std::uint64_t> output_offset_{0};
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
