#pragma once
#include "../support/types.hpp"


namespace nustar::star2711b::solo::detail {

class ReadStatistics {
public:
    std::vector<std::string> names;
    enum {
        noUnmapped,
        noNoFeature,
        MultiFeature,
        subMultiFeatureMultiGenomic,
        noTooManyWLmatches,
        noMMtoWLwithoutExact,
        yesWLmatch,
        yessubWLmatchExact,
        yessubWLmatch_UniqueFeature,
        yesCellBarcodes,
        yesUMIs,
        nStats
    };
    UnsignedWord V[nStats];
    ReadStatistics() {
        names = {"noUnmapped",
                 "noNoFeature",
                 "MultiFeature",
                 "subMultiFeatureMultiGenomic",
                 "noTooManyWLmatches",
                 "noMMtoWLwithoutExact",
                 "yesWLmatch",
                 "yessubWLmatchExact",
                 "yessubWLmatch_UniqueFeature",
                 "yesCellBarcodes",
                 "yesUMIs"};
        for (std::uint32_t ii = 0; ii < nStats; ii++)
            V[ii] = 0;
    };

    UnsignedWord numInvalidBarcodes() {
        return V[noTooManyWLmatches] + V[noMMtoWLwithoutExact];
    };

    UnsignedWord numMappedToTranscriptome() {
        return V[yesWLmatch];
    };

    /*
    void calcUnique(bool multYes)
    {
        if ( multYes  ) {//multi-genic reads were counted in yesWLmatch
            V[yessubWLmatch_UniqueFeature] = V[yesWLmatch] - V[MultiFeature];
        } else {//no multi-genic reads
            V[yessubWLmatch_UniqueFeature] =  V[yesWLmatch];
        };
    };
    */

    UnsignedWord numMappedToTranscriptomeUnique() {
        return V[yessubWLmatch_UniqueFeature];
    };

    double numSequencingSaturation() {
        return 1.0 -
               double(V[yesUMIs]) /
                   V[yessubWLmatch_UniqueFeature]; // yesUMIs is calculated for unqiue-gene reads
    };
};

} // namespace nustar::star2711b::solo::detail
