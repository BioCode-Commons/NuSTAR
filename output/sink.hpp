#pragma once

#include <optional>

#include <chrono>
#include <cstddef>
#include <filesystem>
#include <memory>
#include <span>
#include <vector>
#include "../sam/encode.hpp"
#include "../reference/view.hpp"
#include "../sam/config.hpp"
#include "../statistics/interface.hpp"

namespace nustar::star2711b {

// A movable value owns all bytes until the writer acknowledges completion.
struct OutputBatch {
    BatchId batch;
    std::vector<EncodedFragment> fragments;
};
// A backend-owned temporary chunk. It contains no borrowed worker buffers.
struct StagedOutput {
    explicit StagedOutput(BatchId id) : batch(id) {}
    virtual ~StagedOutput() = default;
    BatchId batch;
};
class IOutputSink {
public:
    virtual ~IOutputSink() = default;
    virtual Result<EncodedFragment> encode(const FragmentId& id, const ProjectedOutput& output) const {
        return encodeSam(id, output);
    }
    // Called after optional Solo counting, when corrected barcode/UMI values exist.
    virtual Status finalize() { return std::monostate{}; }
    virtual Status begin(std::span<const std::byte> sam_header) = 0;
    virtual Status submit(OutputBatch&&) = 0; // Blocking backpressure is permitted.
    // Optional STAR-style chunk spooling: stage is thread-safe, while commit
    // remains serial and ordered. Workers may free all SAM bytes after stage.
    virtual bool supportsStaging() const noexcept { return false; }
    virtual Result<std::unique_ptr<StagedOutput>> stage(OutputBatch&&) {
        return Failure{ErrorKind::Config,Stage::Commit,{},1,"output staging is unsupported"};
    }
    virtual Status submitStaged(std::unique_ptr<StagedOutput>) {
        return Failure{ErrorKind::Config,Stage::Commit,{},1,"output staging is unsupported"};
    }
    virtual Status finish(std::span<const std::byte> sj_table,
        const WorkerStatistics& run_statistics) = 0;
    virtual void cancel() noexcept = 0;
};

struct OutputFiles {
    std::filesystem::path sam;
    std::chrono::system_clock::time_point job_started, mapping_started;
    std::optional<std::string> outFileNamePrefix;
    bool sam_enabled=true;
};
class IOutputSinkBackend {
public:
    virtual ~IOutputSinkBackend() = default;
    virtual Result<std::unique_ptr<IOutputSink>> open(const OutputFiles&) const = 0;
};
const IOutputSinkBackend& bufferedSamSinkBackend();

} // namespace nustar::star2711b
