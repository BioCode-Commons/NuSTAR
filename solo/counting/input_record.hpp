#pragma once

#include "record_types.hpp"
#include "../support/types.hpp"
#include <array>
#include <fstream>
#include <vector>

namespace nustar::star2711b::solo::detail {

bool soloInputFeatureUMI(std::fstream* strIn, std::int32_t featureType, bool readInfoYes,
                         std::array<std::vector<UnsignedWord>, 2>& sjAll, UnsignedWord& iread,
                         std::int32_t& cbmatch, std::uint32_t& feature, UnsignedWord& umi,
                         std::vector<std::uint32_t>& featVecU32, ReadFlags& readFlagCounts);

} // namespace nustar::star2711b::solo::detail
