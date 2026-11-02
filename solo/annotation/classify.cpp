#include "index.hpp"
#include "overlap.hpp"
#include "read_annotations.hpp"
#include "../support/algorithms.hpp"
#include <bitset>

namespace nustar::star2711b::solo::detail {

int alignToTranscript(Alignment& aG, UnsignedWord trS1, std::uint16_t exN1, std::uint32_t* exSE1,
                      std::uint32_t* exLenCum1, std::array<std::uint32_t, 2>& distTrEnds) {
    bool alignIntronic = false;
    bool alignExonic = false;
    bool alignSpansExonIntr = false;
    bool alignSJconcordant = true;

    // we assume that align is fully contained in the transcript, i.e. alignStart>=trStart,
    // alignEnd<=trEnd find exon that overlaps beginning of the read

    // TODO
    // iab=0;
    // distTSS=aG.exons[iab][EX_G]-trS1-exSE1[2*ex1]+exLenCum1[ex1];
    // iab=aG.nExons-1;
    // distTTS=trLen[tr1]-(gthaG.exons[iab][EX_G]-trS1-exSE1[2*ex1]+exLenCum1[ex1]+aG.exons[iab][EX_L]);
    // if trStr1==2: TTS=trLen[tr1]-TTS, TSS=trLen[tr1]-TSS

    for (std::uint32_t iab = 0, ex1 = 0, bS = 0, bE = 0, eE = 0, enS = 0; iab < aG.nExons;
         iab++) { // scan through all blocks of the align

        UnsignedWord bEprev = bE;

        if (aG.exons[iab][EX_G] < trS1) // this can happen for PE reads, when the 2nd mate protrudes
                                        // to the left of the first
            return -1;

        bS = (std::uint32_t)(aG.exons[iab][EX_G] - trS1); // block start
        bE = bS + aG.exons[iab][EX_L] - 1;                // block end

        if (iab == 0 || aG.canonSJ[iab - 1] == -3) { // start of alig, or jump to another mate
            if (!binarySearch_leLeft<std::uint32_t>(bS, exSE1, 2 * exN1, ex1))
                return -1;                     // bS is outside of exons for this transcript
            ex1 = ex1 / 2;                     // ex1 index, with alignStart>=ex1start
        } else if (aG.canonSJ[iab - 1] >= 0) { // splice junction
            if (bEprev == eE && bS == enS) {   // eE and enS are still from the old ex1
                ++ex1;                         // junction agrees
            } else {
                alignSJconcordant = false;
                if (!binarySearch_leLeft<std::uint32_t>(bS, exSE1, 2 * exN1, ex1))
                    return -1; // bS is outside of exons for this transcript
                ex1 = ex1 / 2; // ex1 index, with alignStart>=ex1start
            };
        };
        eE = exSE1[2 * ex1 + 1];
        enS = ex1 + 1 < exN1 ? exSE1[2 * (ex1 + 1)] : 0; // next exon start

        if (bS <= eE) { // block starts in the ex1 exon
            if (bE > eE) {
                alignSpansExonIntr = true;
                // break;//if ex/in span is detected, no need to check anything else - no true :
                // might still have non-concordant junction
            };
            alignExonic = true;

            if (iab == 0) {
                distTrEnds[0] = exLenCum1[ex1] + bS - exSE1[2 * ex1];
            };
            distTrEnds[1] =
                eE - bE +
                (ex1 == (std::uint32_t)exN1 - 1 ? 0
                                                : exSE1[2 * exN1 - 1] - exSE1[2 * exN1 - 2] + 1 +
                                                      exLenCum1[exN1 - 1] - exLenCum1[ex1 + 1]);

        } else {             // block starts in the intron
            if (bE >= enS) { // block overlaps next exon
                alignSpansExonIntr = true;
                // break;//if ex/in span is detected, no need to check anything else
            };
            alignIntronic = true;
        };
    }; // cycle over align blocks

    if (!alignSJconcordant) // if align has a junction, it's always checked for concordance
        return -1; // even for exon/intron aligns, sjs have to be concordant, otherwise align is not
                   // consistent with this transcript model

    if (alignSpansExonIntr) {
        return AlignVsTranscript::ExonIntronSpan; // align spans exon/intron boundary

    } else if (!alignIntronic) {              // align is purely exonic
        return AlignVsTranscript::Concordant; // align is fully exonic and concordant

    } else { // align has introns
        if (alignExonic) {
            return AlignVsTranscript::ExonIntron; // mixed exonic/intronic align, but no span
        } else {
            return AlignVsTranscript::Intron; // purely intronic align
        };
    };
};

int alignToTranscriptMinOverlap(Alignment& aG, UnsignedWord trS1, std::uint32_t* exSE1,
                                std::uint16_t exN1, std::uint32_t minOverlapMinusOne) {
    bool alignIntronic = false;
    bool alignExonic = false;
    bool alignSpansExonIntr = false;
    bool alignSJconcordant = true;

    // we assume that align is fully contained in the transcript, i.e. alignStart>=trStart,
    // alignEnd<=trEnd find exon that overlaps beginning of the read

    for (std::uint32_t iab = 0; iab < aG.nExons; iab++) { // scan through all blocks of the align

        std::uint32_t bS = (std::uint32_t)(aG.exons[iab][EX_G] - trS1); // block start

        std::uint32_t ex1 =
            binarySearch1<std::uint32_t>(bS, exSE1, 2 * exN1) / 2; // alignStart>=ex1start
        if ((std::uint16_t)ex1 == exN1 - 1) {                      // reached last exon.
            // we assume that the end of align is < end of transcript was checked before this
            // function was called ex1 is positive since align is entirely inside this transcript
            alignExonic = true; // can only be exonic
            break;
        };

        while (iab < aG.nExons - 1 && aG.canonSJ[iab] > -3 &&
               aG.canonSJ[iab] < 0) { // indel, expand the block
            ++iab;
        };

        std::uint32_t bE =
            (std::uint32_t)(aG.exons[iab][EX_G] - trS1 + aG.exons[iab][EX_L] - 1); // block end

        if (bE - bS < minOverlapMinusOne) // block is too short
            continue;

        std::uint32_t eE = exSE1[2 * ex1 + 1];  // exon1 end
        std::uint32_t enS = exSE1[2 * ex1 + 2]; // exon2 start
        std::uint32_t enE = exSE1[2 * ex1 + 3]; // exon2 end

        // bS>=eS always
        if (bS + minOverlapMinusOne <= eE) {     // start is certainly in exon1
            if (bE <= eE + minOverlapMinusOne) { // end is in exon1
                alignExonic = true;
            } else { // end spans into intron1
                alignSpansExonIntr = true;
            };

        } else if (bS + minOverlapMinusOne < enS) { // start is in the intron1
            if (bE >= enS + minOverlapMinusOne) {   // end is certainly in exon2 or intron2
                alignSpansExonIntr = true;
            } else if (bE > eE + minOverlapMinusOne) { // end is certainly in intron1
                if (enS - eE >
                    1000000) { // if intron is too large, do not mark this align as intronic. This
                               // is to match velocyto.py logic. TODO: make it a parameter
                    return -1; // no intronic match since the intron is too large. This prevents
                               // large introns from swallowing small genes inside them
                };
                alignIntronic = true;
            }; // otherwise start and end are too close to exon1 end, no call

        } else {                                 // start is too close to the exon2 start
            if (bE > enE + minOverlapMinusOne) { // end is certainly in intron2
                alignSpansExonIntr = true;
            } else if (bE >= enS + minOverlapMinusOne) { // end is certainly in exon2
                alignExonic = true;
            }; // otherwise start and end are too close to exon2 end, no call
        };

        if (aG.sjYes &&
            (alignIntronic || alignSpansExonIntr)) { // spliced align cannot overlap an intron
            alignSJconcordant = false;
            break;
        };
    }; // cycle over align blocks

    if (!alignSJconcordant) // if align has a junction, it's always checked for concordance
        return -1; // even for exon/intron aligns, sjs have to be concordant, otherwise align is not
                   // consistent with this transcript model

    if (alignSpansExonIntr) {
        return AlignVsTranscript::ExonIntronSpan; // align spans exon/intron boundary
    } else if (!alignIntronic) {                  // align is purely exonic
        return AlignVsTranscript::Concordant;     // align is fully exonic and concordant
    } else {                                      // align has introns
        if (alignExonic) {
            return AlignVsTranscript::ExonIntron; // mixed exonic/intronic align, but no span
        } else {
            return AlignVsTranscript::Intron; // purely intronic align
        };
    };
};

void AnnotationIndex::classifyAlign(Alignment** alignG, UnsignedWord nAlignG,
                                    ReadAnnotations& readAnnot) {
    // readAnnot.transcriptConcordant={};
    // readAnnot.trVelocytoType={};
    // array<bool,AlignVsTranscript::N> reAnn={false};

    ReadAnnotFeature& annFeat = readAnnot.annotFeatures[FeatureTypes::Gene];
    // annFeat.fSet={};
    // annFeat.fAlign = {};
    // annFeat.ovType = 0;

    annFeat.fAlign.resize(nAlignG);
    std::uint32_t reGe = (std::uint32_t)-2;  // so that the first gene can be recorded
    std::bitset<velocytoTypeGeneBits> reAnn; // initialized to 0 (false)

    for (UnsignedWord iag = 0; iag < nAlignG; iag++) {

        Alignment& aG = *alignG[iag];

        // binary search through transcript starts
        std::uint32_t tr1 = binarySearch1a<UnsignedWord>(
            aG.exons[0][EX_G], trS.get(),
            nTr); // tr1 has the maximum transcript start such that it is still <= align start
        if (tr1 == (std::uint32_t)-1)
            continue; // this alignment is outside of range of all transcripts

        UnsignedWord aGend = aG.exons[aG.nExons - 1][EX_G] + aG.exons[aG.nExons - 1][EX_L] -
                             1; // TODO: this estimate does work if 2nd mate end is < 1st mate end

        ++tr1;
        do { // cycle back through all the transcripts
            --tr1;
            if (aGend > trE[tr1] ||
                (context_.options.strand >= 0 &&
                 (trStr[tr1] == 1 ? aG.Str : 1 - aG.Str) !=
                     (std::uint32_t)context_.options
                         .strand)) //!(this transcript contains the read, and has correct strand)
                continue;

            std::array<std::uint32_t, 2> distTrEnds;
            int aStatus = alignToTranscript(aG, trS[tr1], trExN[tr1], exSE.get() + 2 * trExI[tr1],
                                            exLenCum.get() + trExI[tr1], distTrEnds);

            if (aStatus == AlignVsTranscript::Concordant) { // align conforms with the transcript

                // debug
                //     cerr << aG.exons[0][EX_L]+distTrEnds[0]+distTrEnds[1] <<" "<<trLen[tr1]
                //     <<endl;

                UnsignedWord distTTS =
                    (UnsignedWord)(trStr[tr1] == 1 ? distTrEnds[1] : distTrEnds[0]);
                readAnnot.transcriptConcordant.push_back({tr1, (std::uint32_t)distTTS});

                annFeat.fSet.insert(trGene[tr1]); // genes for all alignments
                annFeat.fAlign[iag].insert(trGene[tr1]);
            };

            if ((context_.options.featureYes[FeatureTypes::Velocyto] ||
                 context_.options.featureYes[FeatureTypes::VelocytoSimple]) &&
                nAlignG == 1) { // another calculation for velocyto with minOverlapMinusOne=6

                // 6 is the hard code minOverlapMinusOne, to agree with velocyto's MIN_FLANK=5
                aStatus =
                    alignToTranscriptMinOverlap(aG, trS[tr1], exSE.get() + 2 * trExI[tr1], trExN[tr1], 6);

                if (aStatus < 0)
                    continue; // align is not concordant with this transcript because one of the
                              // align junctions does not match transcript junction

                // calculate reAnn
                if (reGe != (std::uint32_t)-1) { // not multi-mapper

                    if (reGe == (std::uint32_t)-2) // first gene
                        reGe = trGene[tr1];

                    if (reGe != trGene[tr1]) {    //
                        reGe = (std::uint32_t)-1; // marks multi-gene align
                    } else {
                        if (aStatus != AlignVsTranscript::ExonIntronSpan) { // no exon/intron span
                            reAnn.set(AlignVsTranscript::ExonIntronSpan,
                                      true); // meaning of this bit is NoExonIntronSpan
                            reAnn.set(aStatus, true);
                        };
                    };
                };

                std::uint8_t reAnn1 = 0; // initialized to 0 (false)
                reAnn1 = 1 << aStatus;
                if (aStatus == AlignVsTranscript::ExonIntronSpan) { // span
                    reAnn1 = reAnn1 |
                             (1 << AlignVsTranscript::Intron); // span is also considered intronic
                    reAnn1 = reAnn1 |
                             (1 << AlignVsTranscript::Concordant); // span is also considered exonic
                };
                readAnnot.trVelocytoType.push_back({tr1, reAnn1});
            };
        } while (trEmax[tr1] >= aGend && tr1 > 0);
    };

    if (annFeat.fSet.size() > 0)
        annFeat.ovType = ReadAnnotFeature::overlapTypes::exonic;
    // VelocytoSimple logic
    readAnnot.geneVelocytoSimple[0] =
        (reGe + 2 == 0 ? (std::uint32_t)-1 : reGe); //-2 marks no gene, convert to -1 which marks
                                                    //either no gene or multigene - no output
    readAnnot.geneVelocytoSimple[1] = reAnn.to_ulong();
};

} // namespace nustar::star2711b::solo::detail
