#pragma once

#include "../common/types.hpp"

namespace nustar::star2711b {

struct ScoreConfig {
    Score scoreGap, scoreGapNoncan, scoreGapGCAG, scoreGapATAC;
    Score scoreDelOpen, scoreDelBase, scoreInsOpen, scoreInsBase;
    Score scoreStitchSJshift, sjdbScore;
    double scoreGenomicLengthLog2scale;
};

} // namespace nustar::star2711b
