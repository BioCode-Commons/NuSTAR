#pragma once

#include "interface.hpp"

namespace nustar::star2711b {

// STAR 2.7.11b short-read window creation and ordered WA distribution.
// The seed source must deliver complete occurrences in physical SA rank order.
Result<std::shared_ptr<const IWindowBuilder>> makeStarWindowBuilder(
    ReferenceLease reference);

} // namespace nustar::star2711b
