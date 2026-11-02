#pragma once

#include "read_recorder.hpp"
#include <array>

namespace nustar::star2711b::solo::detail {

// Native NuSTAR workers provide records and quality histograms directly.
// No ReadAlign/ReadAlignChunk object or STAR mapping engine is instantiated.
struct WorkerState {
    std::unique_ptr<ReadRecorder> read;
    std::array<std::array<UnsignedWord, 256>, 2> qualHist{};
};

} // namespace nustar::star2711b::solo::detail
