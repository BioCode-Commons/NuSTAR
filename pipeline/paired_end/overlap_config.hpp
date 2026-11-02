#pragma once

#include <cstdint>
#include "../../common/types.hpp"

namespace nustar::star2711b {

struct OverlapConfig {
    std::uint64_t peOverlapNbasesMin; // Zero disables overlap merge/remapping.
    double peOverlapMMp;
};

} // namespace nustar::star2711b
