#include "alignment.hpp"

namespace nustar::star2711b::solo::detail {

void Alignment::extractSpliceJunctions(std::vector<std::array<UnsignedWord, 2>>& sjOut,
                                       bool& annotYes) {
    annotYes = true;
    for (UnsignedWord iex = 0; iex < nExons - 1; iex++) { // record all junctions
        if (canonSJ[iex] >= 0) { // only record junctions, not indels or mate gap
            std::array<UnsignedWord, 2> sj;
            sj[0] = exons[iex][EX_G] + exons[iex][EX_L]; // start
            sj[1] = exons[iex + 1][EX_G] - sj[0];        // gap
            sjOut.push_back(sj);
            if (sjAnnot[iex] == 0)
                annotYes = false; // if one of the SJs is unannoated, annotYes=false
        };
    };
};

UnsignedWord Alignment::chrStartLengthExtended() {
    UnsignedWord start1 = cStart - exons[0][EX_R];
    UnsignedWord length1 =
        exons[nExons - 1][EX_G] + Lread - exons[nExons - 1][EX_R] - exons[0][EX_G] + exons[0][EX_R];

    return (start1 << 32) | length1;
};

} // namespace nustar::star2711b::solo::detail
