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
    // Complete batches may reserve disjoint file ranges and write concurrently.
    // This deliberately does not preserve order between input fragments.
    virtual bool concurrentWrites() const noexcept { return false; }
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
