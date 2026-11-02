#pragma once

#include "types.hpp"

namespace nustar::star2711b {

struct StageMarkers {
    std::uint64_t legacy_map_marker; // Retain STAR marker codes and their precedence.
    bool normal_early_stop; // Algorithmic termination is not a process failure.
    std::vector<Diagnostic> diagnostics; // Includes legacy warning-and-truncation events.
};

} // namespace nustar::star2711b
