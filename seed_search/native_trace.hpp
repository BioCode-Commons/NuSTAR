#pragma once

#include "../index/native_index.hpp"

namespace nustar::star2711b {

namespace star_compat {
struct SeedTraceRecord {
    ReadOffset start;
    std::uint64_t length, repeat_count;
    Orientation search_direction;
    SaRange suffixes;
    MateId mate;
};
struct SeedTraceSummary {
    std::uint64_t expanded_hits;
    std::array<std::uint64_t, 2> unique_multi_hits;
    std::uint64_t stored_min_length, longest_unique_length, longest_unique_legacy_index;
    std::uint64_t longest_multi_length, longest_multi_repeats;
    std::uint64_t min_repeats, min_repeats_length, max_repeats, max_repeats_length;
};
struct SeedTrace {
    std::vector<SeedTraceRecord> pieces;
    SeedTraceSummary summary; // Encoded as a versioned BackendTrace, never required by window.
};
} // namespace star_compat

} // namespace nustar::star2711b
