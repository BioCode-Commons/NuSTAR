#pragma once

#include <cstdint>
#include <limits>
#include <span>

namespace nustar::star2711b {

// Reorder SA-IS suffixes that STAR 2.7.11b treats as equal after chromosome/SJ
// padding or genomeSuffixLengthMax. The coordinate map reflects the current
// four-segment NuSTAR physical sequence layout.
void normalizeStarSuffixRanks(std::span<std::uint32_t> suffixes,
    std::span<const std::uint8_t> sequence, std::uint64_t real_length,
    std::uint64_t synthetic_length,
    std::uint32_t suffix_length_max = std::numeric_limits<std::uint32_t>::max());

} // namespace nustar::star2711b
