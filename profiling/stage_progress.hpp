#pragma once

#include <array>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>

namespace nustar::star2711b {

// Cumulative, per-read stage timing. The current trace tools are single-threaded;
// a future multi-worker runner should give each worker a recorder and merge rows.
enum class TimedStage : std::uint8_t { Seed = 0, Cluster = 1, Stitch = 2 };
constexpr std::uint8_t timedStageBit(TimedStage stage) {
    return static_cast<std::uint8_t>(1U << static_cast<unsigned>(stage));
}

class StageProgress {
public:
    // Zero disables timing; nonzero intervals below 10,000 are rounded up.
    StageProgress(std::filesystem::path path, std::uint64_t every_reads,
                  std::uint8_t available_stages);
    ~StageProgress();
    StageProgress(const StageProgress&) = delete;
    StageProgress& operator=(const StageProgress&) = delete;

    static StageProgress fromEnvironment(std::filesystem::path default_path,
                                         std::uint8_t available_stages);
    bool enabled() const noexcept { return every_reads_ != 0; }
    void beginRead();
    void beginStage(TimedStage stage);
    void endStage(TimedStage stage);
    void endRead();
    void finish();

private:
    using Clock = std::chrono::steady_clock;
    void writeRow(const char* status);
    static std::uint64_t nanoseconds(Clock::duration duration);

    std::ofstream output_;
    std::uint64_t every_reads_ = 0, completed_reads_ = 0;
    std::uint8_t available_stages_ = 0;
    bool active_read_ = false, active_stage_ = false, finished_ = false;
    TimedStage current_stage_ = TimedStage::Seed;
    Clock::time_point wall_start_{}, read_start_{}, stage_start_{};
    std::array<std::uint64_t, 3> read_stage_ns_{}, total_stage_ns_{};
    std::uint64_t total_read_ns_ = 0;
};

} // namespace nustar::star2711b
