#include "sequence.hpp"
#include <assert.h>

namespace nustar::star2711b::solo::detail {

char convertNt01234(const char R0) {//transform sequence  from ACGT into 0-1-2-3 code
    switch(R0)
    {
        case('a'):
        case('A'):
            return 0;
            break;
        case('c'):
        case('C'):
            return 1;
            break;
        case('g'):
        case('G'):
            return 2;
            break;
        case('t'):
        case('T'):
            return 3;
            break;
        default:
            return 4;
    };
};


SignedWord convertNuclStrToInt64(const std::string S, UnsignedWord& intOut) {
    intOut = 0;
    SignedWord posN = -1;
    for (UnsignedWord ii = 0; ii < S.size(); ii++) {
        UnsignedWord nt = (UnsignedWord)convertNt01234(S[ii]);
        if (nt > 3) { // N
            if (posN >= 0)
                return -2; // two Ns
            posN = ii;
            nt = 0;
        };
        intOut = intOut << 2;
        intOut += nt;
        // intOut += nt<<(2*ii);
    };
    return posN;
};

std::string convertNuclInt64toString(UnsignedWord nuclNum, const std::uint32_t L) {
    std::string nuclOut(L, 'N');
    std::string nuclChar = "ACGT";

    for (UnsignedWord ii = 1; ii <= L; ii++) {
        nuclOut[L - ii] = nuclChar[nuclNum & 3];
        nuclNum = nuclNum >> 2;
    };

    return nuclOut;
};

std::uint32_t localAlignHammingDist(const std::string& text, const std::string& query,
                                    std::uint32_t& pos) {
    std::uint32_t distBest = query.size();
    if (text.size() < query.size()) { // query is longer than text, no match
        return text.size() + 1;
    };
    for (std::uint32_t ii = 0; ii < text.size() - query.size() + 1; ii++) {
        std::uint32_t dist1 = 0;
        for (std::uint32_t jj = 0; jj < query.size(); jj++) {
            if (query[jj] != 'N' &&
                text[jj + ii] != query[jj]) { // N in query does not count as mismatch
                ++dist1;
            };
        };
        if (dist1 < distBest) {
            distBest = dist1;
            pos = ii;
        };
    };
    return distBest;
};
} // namespace nustar::star2711b::solo::detail
