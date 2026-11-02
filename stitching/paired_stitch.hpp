#pragma once

#include "paired_interface.hpp"

namespace nustar::star2711b {

Result<std::shared_ptr<const IPairedStitcher>> makeStarPairedStitcher(
    ReferenceLease reference);

} // namespace nustar::star2711b
