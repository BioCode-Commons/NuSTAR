#pragma once

#include "../support/types.hpp"

namespace nustar::star2711b::solo::detail {

/**  general basic functions to operate on sequences, no classes
 */

std::uint32_t localAlignHammingDist(const std::string& text, const std::string& query,
                                    std::uint32_t& pos);

SignedWord convertNuclStrToInt64(const std::string S, UnsignedWord& intOut);
std::string convertNuclInt64toString(const UnsignedWord nuclNum, const std::uint32_t L);

} // namespace nustar::star2711b::solo::detail
