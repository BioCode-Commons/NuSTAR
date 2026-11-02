#pragma once

#include "alignment.hpp"
#include "overlap.hpp"
#include "read_annotations.hpp"
#include "../configuration/context.hpp"
#include "../support/types.hpp"
#include <memory>
#include <set>

namespace nustar::star2711b::solo::detail {

class AnnotationIndex {
public:
    std::string trInfoDir;

    std::vector<std::string> trID, geID, geName, geBiotype; // transcript/gene IDs
    std::uint32_t nTr, nGe;                                 // number of transcript/genes

    std::unique_ptr<UnsignedWord[]> trS, trE, trEmax; // transcripts start,end,end-max

    std::uint32_t nEx;                      // number of exons
    std::unique_ptr<std::uint16_t[]> trExN; // number of exons per transcript
    std::unique_ptr<std::uint32_t[]> trExI; // index of the first exon for each transcript in exSE
    std::unique_ptr<std::uint8_t[]> trStr;  // transcript strand
    std::unique_ptr<std::uint32_t[]> exSE;  // exons start/end
    std::unique_ptr<std::uint32_t[]> exLenCum; // cumulative length of previous exons
    std::unique_ptr<std::uint32_t[]> trGene;   // transcript to gene correspondence
    std::unique_ptr<std::uint32_t[]> trLen;    // transcript lengths

    struct {                                        // exon-gene structure for GeneCounts
        UnsignedWord nEx;                           // number of exons/genes
        std::unique_ptr<UnsignedWord[]> s, e, eMax; // exon start/end
        std::unique_ptr<std::uint8_t[]> str;        // strand
        std::unique_ptr<std::uint32_t[]> g, t;      // gene/transcript IDs
    } exG;

    struct { // geneFull structure
        std::unique_ptr<UnsignedWord[]> s, e, eMax;
        std::unique_ptr<std::uint8_t[]> str;
        std::unique_ptr<std::uint32_t[]> g;
    } geneFull;

    // methods:
    ~AnnotationIndex() = default;
    AnnotationIndex(const AnnotationIndex&) = delete;
    AnnotationIndex& operator=(const AnnotationIndex&) = delete;
    AnnotationIndex(
        Context& context_in); // create transcriptome structure, load and initialize parameters
                                                // in RA to transcriptomic in RAtr
    void geneFullAlignOverlap(UnsignedWord nA, Alignment** aAll, std::int32_t strandType,
                              ReadAnnotFeature& annFeat);
    void geneFullAlignOverlap_ExonOverIntron(UnsignedWord nA, Alignment** aAll,
                                             std::int32_t strandType, ReadAnnotFeature& annFeat,
                                             ReadAnnotFeature& annFeatGeneConcordant);
    void classifyAlign(Alignment** alignG, UnsignedWord nAlignG, ReadAnnotations& readAnnot);
    void alignExonOverlap(UnsignedWord nA, Alignment** aAll, std::int32_t strandType,
                          ReadAnnotFeature& annFeat);

private:
    Context& context_; // normal "genomic" parameters
};

} // namespace nustar::star2711b::solo::detail
