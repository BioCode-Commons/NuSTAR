#pragma once

#include <cstdint>
#include <vector>
#include "../reference/genome.hpp"

namespace nustar::star2711b {

// Four physical segments: forward real, reverse-complement real, forward SJ,
// reverse-complement SJ, each followed by a 20-base padding gap.
// shift=1 reserves zero for SAIS's sentinel; shift=0 uses STAR base codes.
std::vector<std::uint8_t> makeStarPhysicalSequence(
    const GenomeReferenceData& reference, std::uint8_t shift);

} // namespace nustar::star2711b
