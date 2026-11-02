#include "index.hpp"
#include "read_annotations.hpp"
#include "../support/algorithms.hpp"

namespace nustar::star2711b::solo::detail {

void AnnotationIndex::geneFullAlignOverlap(UnsignedWord nA, Alignment** aAll,
                                           std::int32_t strandType, ReadAnnotFeature& annFeat) {
    // annFeat.fSet={};
    // annFeat.fAlign = {};
    // annFeat.ovType = 0; //exonic/intronic determination is not done

    annFeat.fAlign.resize(nA);
    for (std::uint32_t iA = 0; iA < nA; iA++) {
        Alignment& a = *aAll[iA]; // one unique alignment only

        SignedWord gi1 = -1;

        for (SignedWord ib = a.nExons - 1; ib >= 0;
             ib--) { // scan through all blocks of the alignments

            UnsignedWord be1 = a.exons[ib][EX_G] + a.exons[ib][EX_L] - 1; // end of the block
            gi1 = binarySearch1a<UnsignedWord>(be1, geneFull.s.get(), (std::int32_t)nGe);

            while (gi1 >= 0 &&
                   geneFull.eMax[gi1] >= a.exons[ib][EX_G]) { // these exons may overlap this block
                if (geneFull.e[gi1] >= a.exons[ib][EX_G]) {   // this gene overlaps the block
                    std::int32_t str1 = geneFull.str[gi1] == 1 ? a.Str : 1 - a.Str;
                    if (strandType == -1 || strandType == str1) {
                        annFeat.fSet.insert(geneFull.g[gi1]);
                        annFeat.fAlign[iA].insert(geneFull.g[gi1]);
                    };
                };
                --gi1; // go to the previous gene
            };
        };
    };

    
};

} // namespace nustar::star2711b::solo::detail
