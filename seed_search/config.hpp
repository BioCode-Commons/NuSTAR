#pragma once

#include <cstdint>
#include "../common/types.hpp"

namespace nustar::star2711b {

struct SeedConfig {
    std::uint64_t seedSearchStartLmax, seedSearchLmax;
    double seedSearchStartLmaxOverLread;
    std::uint64_t seedMultimapNmax, seedPerReadNmax;
    std::uint64_t seedSplitMin, seedMapMin;
    std::uint64_t maxNsplit; // Legacy internal limit; minLmap appears only in stale comments.
    std::uint64_t outFilterMatchNmin; // Used by seed-stage early termination.
};

} // namespace nustar::star2711b
