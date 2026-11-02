#pragma once

#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <vector>
#include "../sam/model.hpp"
#include "../stitching/candidate.hpp"
#include "../selection/model.hpp"
#include "../reads/model.hpp"

namespace nustar::star2711b {

struct StatsDelta {
    std::uint64_t reads, read_bases, mapped_unique, mapped_multiple;
    std::uint64_t mapped_bases, mismatches, insertion_count, deletion_count;
    std::uint64_t insertion_bases, deletion_bases;
    double mapped_portion;
    std::array<std::uint64_t, 7> splices_by_motif;
    std::uint64_t annotated_splices;
    std::uint64_t unmapped_other, unmapped_short, unmapped_mismatch, unmapped_multi, unmapped_all;
};
struct WorkerStatistics {
    StatsDelta totals;
    std::array<std::array<std::uint64_t, 256>, 2> quality_histogram;
};
struct QualityCount { std::uint8_t value; std::uint32_t count; };
struct QualityHistogramDelta {
    std::array<std::vector<QualityCount>, 2> by_mate; // Nonzero bins, ascending byte value.
};
Result<StatsDelta> deriveStatistics(const MappingRead&, const CandidateSet&,
    const SelectionResult&, const FilterDecision&, const OutputDisposition&);
void addStatistics(StatsDelta& total, const StatsDelta& delta) noexcept;
QualityHistogramDelta deriveQualityHistogram(
    std::span<const std::string_view> quality_by_mate);
std::string formatStarFinalStatistics(const StatsDelta&,
    std::chrono::system_clock::time_point started,
    std::chrono::system_clock::time_point mapping_started,
    std::chrono::system_clock::time_point finished);

} // namespace nustar::star2711b
