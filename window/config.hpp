#pragma once

#include <cstdint>
#include "../common/types.hpp"

namespace nustar::star2711b {

struct WindowConfig {
    std::uint64_t winAnchorMultimapNmax;
    std::uint32_t winBinNbits, winBinChrNbits;
    std::uint64_t winBinN, winAnchorDistNbins, winFlankNbins;
    std::uint64_t seedPerWindowNmax, seedNoneLociPerWindow;
    std::uint64_t alignWindowsPerReadNmax;
    double winReadCoverageRelativeMin;
    std::uint64_t winReadCoverageBasesMin;
};

} // namespace nustar::star2711b
