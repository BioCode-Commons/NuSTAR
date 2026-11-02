#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>
#include "../common/types.hpp"

namespace nustar::star2711b {

// Config field names intentionally retain STAR parameter spelling.
// Objects contain resolved values, not an alternative set of hard-coded defaults.
struct ClipEndConfig {
    std::uint64_t fixed_bases;
    std::string adapter;
    double adapter_mismatch_fraction;
    std::uint64_t bases_after_adapter;
};
struct PrepareCommonConfig {
    std::string clipAdapterType; // Validated STAR token, including Hamming/CellRanger4.
    std::vector<char> readNameSeparator;
    std::uint32_t readQualityScoreBase;
    int outQSconversionAdd; // FASTQ quality offset, clamped to ASCII [33,126] as in readLoad.
};
struct SinglePrepareConfig {
    PrepareCommonConfig common;
    std::array<ClipEndConfig, 2> clip; // [5-prime, 3-prime].
};
struct PairedPrepareConfig {
    PrepareCommonConfig common;
    std::array<std::array<ClipEndConfig, 2>, 2> clip; // [mate][5-prime, 3-prime].
};

} // namespace nustar::star2711b
