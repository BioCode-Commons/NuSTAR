#pragma once

#include "interface.hpp"
#include "../reference/view.hpp"

namespace nustar::star2711b {

Result<std::vector<std::byte>> encodeSjTable(const FinalJunctions&,
    const ReferenceMetadata&); // Exactly nine columns, no header, LF line endings.

} // namespace nustar::star2711b
