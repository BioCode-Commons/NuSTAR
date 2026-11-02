#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include "../common/types.hpp"

namespace nustar::star2711b {

struct SelectionConfig {
    Score outFilterMultimapScoreRange;
    std::uint64_t outFilterMultimapNmax;
    std::optional<std::uint64_t> outSAMmultNmax; // Null means the legacy unlimited sentinel.
    std::string outMultimapperOrder, outSAMprimaryFlag;
};
struct ReadFilterConfig {
    Score outFilterScoreMin;
    double outFilterScoreMinOverLread;
    std::uint64_t outFilterMatchNmin, outFilterMismatchNmax;
    double outFilterMatchNminOverLread;
    double outFilterMismatchNoverLmax, outFilterMismatchNoverReadLmax;
};

} // namespace nustar::star2711b
