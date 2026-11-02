#pragma once

#include <array>
#include <cstdint>
#include <string>
#include "../common/types.hpp"
#include "score.hpp"

namespace nustar::star2711b {

struct StitchCommonConfig {
    std::uint64_t alignIntronMin, alignIntronMax;
    std::uint64_t alignSJoverhangMin, alignSJDBoverhangMin;
    std::array<std::int32_t, 4> alignSJstitchMismatchNmax;
    std::uint64_t alignSplicedMateMapLmin;
    double alignSplicedMateMapLminOverLmate;
    std::uint64_t alignTranscriptsPerWindowNmax, alignTranscriptsPerReadNmax;
    std::string alignEndsType;
    bool alignSoftClipAtReferenceEnds;
    bool insertion_flush_right;
    std::string outFilterIntronMotifs, outFilterIntronStrands, outSAMstrandField;
    double outFilterMismatchNoverLmax;
    Score outFilterMultimapScoreRange; // Also used while retaining window candidates.
};
struct SingleStitchConfig {
    StitchCommonConfig common;
    std::array<bool, 2> force_end_extension; // [5-prime, 3-prime].
};
struct PairedStitchConfig {
    StitchCommonConfig common;
    std::array<std::array<bool, 2>, 2> force_end_extension;
    std::uint64_t alignMatesGapMax, protrusion_bases_max;
    bool protrusion_concordant_pair;
};

} // namespace nustar::star2711b
