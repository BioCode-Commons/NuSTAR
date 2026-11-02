#pragma once

#include <cstdint>
#include "types.hpp"

namespace nustar::star2711b {

struct CompatibilityLimits {
    std::uint64_t max_combined_read_length; // DEF_readSeqLengthMax, not per-mate length.
    std::uint64_t max_name_length, max_blocks, max_multimaps;
    std::uint32_t max_sj_repeat_search;
};

} // namespace nustar::star2711b
