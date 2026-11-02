#pragma once

#include "single_interface.hpp"

namespace nustar::star2711b {

// STAR 2.7.11b short-read, single-mate window stitching. The caller supplies
// windows in STAR WA order; window creation and seed distribution are separate.
Result<std::shared_ptr<const ISingleStitcher>> makeStarSingleStitcher(
    ReferenceLease reference);

} // namespace nustar::star2711b
