#pragma once

#include "view.hpp"

namespace nustar::star2711b {

std::vector<BaseCode> makeSyntheticJunctionSequence(
    std::span<const SjdbEntry> junctions, std::span<const BaseCode> forward,
    std::uint64_t overhang);

} // namespace nustar::star2711b
