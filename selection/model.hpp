#pragma once

#include <cstdint>
#include <optional>
#include <random>
#include <string>
#include <vector>
#include "../stitching/candidate.hpp"

namespace nustar::star2711b {

struct SelectedAlignment {
    CandidateId candidate;
    bool primary;
    ReadOffset original_read_start;
    std::uint64_t chromosome_relative_start;
};
struct SelectionResult {
    std::vector<SelectedAlignment> ordered; // Includes candidates later excluded by SAM count limits.
    std::optional<CandidateId> best;
    Score best_score;
    std::uint64_t candidate_count;
    bool normalized_coordinates_valid; // False for early excessive-multimap return.
};
enum class MappingClass { Unique, Multiple, Unmapped };
enum class UnmappedReason { None, NoGoodWindow, TooShort, TooManyMismatches, TooManyLoci };
struct FilterDecision {
    MappingClass mapping;
    UnmappedReason reason;
    int legacy_unmap_type; // -1/-2 mapped; 0..3 rejected; output may separately produce 4.
};
struct LegacyRandomState {
    std::mt19937 engine;
    std::uniform_real_distribution<double> uniform_0_1;
};
struct RandomSnapshot {
    std::string engine_state, distribution_state, standard_library;
};
RandomSnapshot snapshotRandom(const LegacyRandomState&);
Status restoreRandom(const RandomSnapshot&, LegacyRandomState&);

} // namespace nustar::star2711b
