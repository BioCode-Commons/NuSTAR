#pragma once
#include "../support/types.hpp"

namespace nustar::star2711b::solo::detail {

// Only the selected candidate values consumed by gene/SJ/SmartSeq counting.
// The mapper owns scores, variants, CIGAR and clipping. This is not a mapper.
struct Alignment {
    // Reused scratch: Worker fills every active block and boundary before use.
    // Avoid zeroing the full fixed-capacity arrays when growing the workspace.
    Alignment() noexcept {}
    UnsignedWord exons[MAX_N_EXONS][EX_SIZE];
    int canonSJ[MAX_N_EXONS];
    std::uint8_t sjAnnot[MAX_N_EXONS];
    UnsignedWord nExons;
    UnsignedWord Chr, Str;
    UnsignedWord cStart, Lread, readNmates;
    bool sjYes;

    void extractSpliceJunctions(std::vector<std::array<UnsignedWord, 2>>& junctions,
                               bool& all_annotated);
    UnsignedWord chrStartLengthExtended();
};

} // namespace nustar::star2711b::solo::detail
