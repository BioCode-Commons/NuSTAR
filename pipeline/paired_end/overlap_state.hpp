#pragma once

#include <array>
#include <cstdint>
#include "../../stitching/candidate.hpp"

namespace nustar::star2711b {

struct OverlapState {
    bool attempted, accepted;
    std::uint64_t overlapping_bases;
    std::array<std::uint64_t, 2> mate_starts_in_merged;
};

} // namespace nustar::star2711b
