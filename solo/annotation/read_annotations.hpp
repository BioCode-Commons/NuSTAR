#pragma once
#include <set>

#include "feature_types.hpp"
#include "../counting/record_types.hpp"
#include "../support/types.hpp"

namespace nustar::star2711b::solo::detail {

class ReadAnnotFeature { // annotations for one feature
public:
    std::set<std::uint32_t> fSet;                // set of genes for this read
    std::vector<std::set<std::uint32_t>> fAlign; // gene for each alignment of this read
    std::uint32_t ovType;
    enum overlapTypes {
        none,
        exonic,
        exonicAS,
        exonic50p,
        exonic50pAS,
        intronic,
        intronicAS,
        intergenic,
        N
    };

    void reset() {
        fSet = {};
        fAlign = {};
        ovType = 0;
    };
};

class ReadAnnotations {
public:
    // set<uint32> geneFull, geneFull_Ex50pAS, geneFull_ExonOverIntron, geneConcordant;

    std::array<ReadAnnotFeature, FeatureTypes::N> annotFeatures;

    std::vector<std::array<std::uint32_t, 2>> transcriptConcordant;
    std::vector<std::int32_t> geneExonOverlap;
    std::array<std::uint32_t, 2>
        geneVelocytoSimple; // first element is gene, then counts of transcript types
    std::vector<TranscriptType>
        trVelocytoType; // first element is gene, then counts of transcript types

    void reset() {
        transcriptConcordant = {};
        trVelocytoType = {};
        for (auto& raf : annotFeatures)
            raf.reset();
    };
};

} // namespace nustar::star2711b::solo::detail
