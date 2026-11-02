#include "stage_progress.hpp"

#include <cerrno>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>

namespace nustar::star2711b {
namespace {

constexpr std::uint64_t kMinimumProgressReads = 10000;

std::uint64_t intervalFromEnvironment() {
    const char* value = std::getenv("NUSTAR_STAGE_TIMING_EVERY");
    if (!value || !*value) return 0;
    if (*value == '-') throw std::invalid_argument("NUSTAR_STAGE_TIMING_EVERY must be nonnegative");
    errno = 0;
    char* end = nullptr;
    const auto parsed = std::strtoull(value, &end, 10);
    if (errno == ERANGE || !end || *end != '\0' ||
        parsed > std::numeric_limits<std::uint64_t>::max())
        throw std::invalid_argument("invalid NUSTAR_STAGE_TIMING_EVERY");
    return static_cast<std::uint64_t>(parsed);
}

} // namespace

StageProgress::StageProgress(std::filesystem::path path, std::uint64_t every_reads,
                             std::uint8_t available_stages)
    : every_reads_(every_reads == 0 ? 0 :
                   (every_reads < kMinimumProgressReads ?
                    kMinimumProgressReads : every_reads)),
      available_stages_(available_stages) {
    if (!enabled()) return;
    output_.open(path);
    if (!output_) throw std::runtime_error("cannot open stage timing file: " + path.string());
    output_ << "status\treads\telapsed_wall_s\tseed_s\tcluster_s\tstitch_s"
               "\tcluster_stitch_s\tread_total_s\tother_s\n";
    output_ << std::fixed << std::setprecision(6);
    wall_start_ = Clock::now();
}

StageProgress::~StageProgress() {
    try { finish(); } catch (...) {}
}

StageProgress StageProgress::fromEnvironment(std::filesystem::path default_path,
                                             std::uint8_t available_stages) {
    const auto every = intervalFromEnvironment();
    if (every != 0) {
        if (const char* path = std::getenv("NUSTAR_STAGE_TIMING_FILE"); path && *path)
            default_path = path;
    }
    return StageProgress(std::move(default_path), every, available_stages);
}

std::uint64_t StageProgress::nanoseconds(Clock::duration duration) {
    return static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count());
}

void StageProgress::beginRead() {
    if (!enabled()) return;
    if (active_read_) throw std::logic_error("timing read already active");
    active_read_ = true;
    read_stage_ns_ = {};
    read_start_ = Clock::now();
}

void StageProgress::beginStage(TimedStage stage) {
    if (!enabled()) return;
    if (!active_read_ || active_stage_ || !(available_stages_ & timedStageBit(stage)))
        throw std::logic_error("invalid timing stage start");
    current_stage_ = stage;
    active_stage_ = true;
    stage_start_ = Clock::now();
}

void StageProgress::endStage(TimedStage stage) {
    if (!enabled()) return;
    if (!active_stage_ || current_stage_ != stage)
        throw std::logic_error("invalid timing stage end");
    const auto now = Clock::now();
    read_stage_ns_[static_cast<unsigned>(stage)] += nanoseconds(now-stage_start_);
    active_stage_ = false;
}

void StageProgress::endRead() {
    if (!enabled()) return;
    if (!active_read_ || active_stage_) throw std::logic_error("invalid timing read end");
    const auto now = Clock::now();
    total_read_ns_ += nanoseconds(now-read_start_);
    for (std::size_t i = 0; i < 3; ++i) total_stage_ns_[i] += read_stage_ns_[i];
    active_read_ = false;
    ++completed_reads_;
    if (completed_reads_ % every_reads_ == 0) writeRow("progress");
}

void StageProgress::writeRow(const char* status) {
    const auto seconds = [](std::uint64_t ns) { return double(ns) / 1e9; };
    output_ << status << '\t' << completed_reads_ << '\t'
            << seconds(nanoseconds(Clock::now()-wall_start_));
    for (unsigned i = 0; i < 3; ++i) {
        output_ << '\t';
        if (available_stages_ & (1U << i)) output_ << seconds(total_stage_ns_[i]);
        else output_ << "NA";
    }
    const auto cluster_and_stitch = timedStageBit(TimedStage::Cluster) |
                                    timedStageBit(TimedStage::Stitch);
    output_ << '\t';
    if ((available_stages_ & cluster_and_stitch) == cluster_and_stitch)
        output_ << seconds(total_stage_ns_[1] + total_stage_ns_[2]);
    else output_ << "NA";
    output_ << '\t' << seconds(total_read_ns_) << '\t';
    if ((available_stages_ & 7U) == 7U) {
        const auto timed = total_stage_ns_[0]+total_stage_ns_[1]+total_stage_ns_[2];
        output_ << seconds(total_read_ns_ >= timed ? total_read_ns_-timed : 0);
    } else output_ << "NA";
    output_ << '\n';
    output_.flush();
    if (!output_) throw std::runtime_error("cannot write stage timing file");
}

void StageProgress::finish() {
    if (!enabled() || finished_) return;
    if (active_read_) throw std::logic_error("timing read still active");
    writeRow("final");
    finished_ = true;
}

} // namespace nustar::star2711b
