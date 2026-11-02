#pragma once
#include "../support/types.hpp"


namespace nustar::star2711b::solo::detail {

class BarcodeStatistics {
public:
    std::vector<std::string> names;
    enum {
        noNoAdapter,
        noNoUMI,
        noNoCB,
        noNinCB,
        noNinUMI,
        noUMIhomopolymer,
        noNoWLmatch,
        noTooManyMM,
        noTooManyWLmatches,
        yesWLmatchExact,
        yesOneWLmatchWithMM,
        yesMultWLmatchWithMM,
        nStats
    };
    UnsignedWord V[nStats];
    BarcodeStatistics() {
        names = {"noNoAdapter",
                 "noNoUMI",
                 "noNoCB",
                 "noNinCB",
                 "noNinUMI",
                 "noUMIhomopolymer",
                 "noNoWLmatch",
                 "noTooManyMM",
                 "noTooManyWLmatches",
                 "yesWLmatchExact",
                 "yesOneWLmatchWithMM",
                 "yesMultWLmatchWithMM"};
        for (std::uint32_t ii = 0; ii < nStats; ii++)
            V[ii] = 0;
    };

    UnsignedWord numInvalidBarcodes() {
        UnsignedWord n = 0;
        for (std::uint32_t ii = 0; ii < 9; ii++) // first eight number record all invalid barcodes
            n += V[ii];

        return n;
    };
};

} // namespace nustar::star2711b::solo::detail
