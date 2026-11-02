#pragma once

#include <array>
#include <cstdint>
#include <vector>
#include "../common/types.hpp"

namespace nustar::star2711b {

struct SjOutputConfig {
    bool enabled;
    bool unique_reads_only; // outSJfilterReads == Unique.
    std::array<std::int32_t, 4> outSJfilterOverhangMin;
    std::array<std::int32_t, 4> outSJfilterCountUniqueMin;
    std::array<std::int32_t, 4> outSJfilterCountTotalMin;
    std::array<std::int32_t, 4> outSJfilterDistToOtherSJmin;
    std::vector<std::int32_t> outSJfilterIntronMaxVsReadN;
    std::uint64_t limitOutSJoneRead, limitOutSJcollapsed;
};

} // namespace nustar::star2711b
