#include "../barcode/sequence.hpp"
#include "feature_counter.hpp"
#include "record_types.hpp"
#include "../support/memory.hpp"
#include "../support/streams.hpp"
#include "../support/time.hpp"
#include <bitset>
#include <unordered_map>

namespace nustar::star2711b::solo::detail {

void FeatureCounter::countVelocyto() { // velocyto counting gets info from Gene counting
    time_t rawTime;

    nReadPerCB.resize(nCB);

    std::vector<std::unordered_map<UmiCode, std::vector<TranscriptType>>> cuTrTypes(nCB);
    for (std::uint32_t ii = 0; ii < nCB; ii++)
        cuTrTypes[ii].reserve(readFeatSum->cbReadCount[ii] > 100
                                  ? readFeatSum->cbReadCount[ii]
                                  : readFeatSum->cbReadCount[ii] / 5); // heuristics...

    time(&rawTime);
    context_.inOut->logMain << timeMonthDayTime(rawTime)
                            << " ... Velocyto counting: allocated arrays" << std::endl;
    for (int iThread = 0; iThread < context_.runThreadN;
         iThread++) { // TODO: this can be parallelized
        std::fstream* streamReads = readFeatAll[iThread]->streamReads;
        streamReads->flush();
        streamReads->seekg(0, std::ios::beg);

        UnsignedWord iread;
        while (*streamReads >> iread) { // until the end of file
            BarcodeCode cb =
                features_all_[options.featureInd[FeatureTypes::Gene]]->readInfo[iread].cb;
            UmiCode umi =
                features_all_[options.featureInd[FeatureTypes::Gene]]->readInfo[iread].umi;
            if (cb == (BarcodeCode)-1 ||
                umi == (UmiCode)-1) { // CB and/or UMI undefined.   TODO: put a filter on CBs here,
                                      // e.g. UMI threshold
                streamReads->ignore((std::uint32_t)-1, '\n');
                continue;
            };

            std::uint32_t iCB = indCBwl[cb];
            nReadPerCB[iCB]++; // simple estimate

            if (cuTrTypes[iCB].count(umi) > 0 &&
                cuTrTypes[iCB][umi]
                    .empty()) { // intersection is empty, no need to load this transcript
                streamReads->ignore((std::uint32_t)-1, '\n');
                continue;
            };

            std::uint32_t nTr;
            *streamReads >> nTr;
            std::vector<TranscriptType> trT(nTr);
            for (auto& tt : trT) {
                std::uint32_t ty;
                *streamReads >> tt.tr >> ty;
                tt.type = (std::uint8_t)ty;
            };

            if (cuTrTypes[iCB].count(umi) == 0) { // 1st entry for this umi
                cuTrTypes[iCB][umi] = trT;
                continue;
            };

            std::uint32_t inew = 0;
            std::vector<TranscriptType> trT1;
            trT1.reserve(cuTrTypes[iCB][umi].size());

            for (std::uint32_t iold = 0; iold < cuTrTypes[iCB][umi].size();
                 iold++) { // intersection of old with new
                while (inew < trT.size() &&
                       cuTrTypes[iCB][umi][iold].tr > trT[inew].tr) // move through the sorted lists
                    ++inew;

                if (inew == trT.size()) // end of trT reached
                    break;

                if (cuTrTypes[iCB][umi][iold].tr == trT[inew].tr) { //
                    trT1.push_back({trT[inew].tr, (std::uint8_t)(cuTrTypes[iCB][umi][iold].type |
                                                                 trT[inew].type)});
                };
            };
            cuTrTypes[iCB][umi] = trT1; // replace with intersection
        };
    };

    time(&rawTime);
    context_.inOut->logMain << timeMonthDayTime(rawTime) << " ... Velocyto counting: finished input"
                            << std::endl;
    nUMIperCB.resize(nCB, 0);
    nGenePerCB.resize(nCB, 0);

    countMatStride = 4; // spliced/unspliced/ambiguous
    countCellGeneUMI.resize(nReadsMapped * countMatStride / 5 +
                            16); // 5 is heuristic, will be resized if needed
    countCellGeneUMIindex.resize(nCB + 1);
    countCellGeneUMIindex[0] = 0;

    for (std::uint32_t iCB = 0; iCB < nCB; iCB++) { // main collapse cycle
        std::map<std::uint32_t, std::array<std::uint32_t, 3>> geneC;
        for (auto& umi : cuTrTypes[iCB]) { // cycle over UMIs
            if (umi.second.empty())        // no transcripts in the intesect
                continue;
            std::uint32_t geneI = annotation_.trGene[umi.second[0].tr];

            // for each transcript, are all UMIs:
            bool exonModel = false;   // purely exonic
            bool intronModel = false; // purely intronic
            bool spanModel = true;    // spanning model for all UMIs
            bool mixedModel = false;  // both intronic and exonic, but not spanning

            for (auto& tt : umi.second) {                 // cycle over transcripts in this UMI
                if (annotation_.trGene[tt.tr] != geneI) { // multigene
                    geneI = (std::uint32_t)-1;
                    break;
                };

                std::bitset<velocytoTypeGeneBits> gV(tt.type);

                mixedModel |=
                    ((gV.test(AlignVsTranscript::Intron) &&
                      gV.test(AlignVsTranscript::Concordant)) ||
                     gV.test(AlignVsTranscript::ExonIntron)) &&
                    !gV.test(AlignVsTranscript::ExonIntronSpan); // has exon and intron, but no span
                spanModel &= gV.test(AlignVsTranscript::ExonIntronSpan);
                exonModel |= gV.test(AlignVsTranscript::Concordant) &&
                             !gV.test(AlignVsTranscript::Intron) &&
                             !gV.test(AlignVsTranscript::ExonIntron); // has only-exons model
                intronModel |= gV.test(AlignVsTranscript::Intron) &&
                               !gV.test(AlignVsTranscript::ExonIntron) &&
                               !gV.test(AlignVsTranscript::Concordant); // has only-introns model,
                                                                        // this includes span
            };

            if (geneI + 1 == 0) // multigene
                continue;

            if (exonModel && !intronModel && !mixedModel) { // all models are only-exons
                geneC[geneI][0]++;                          // spliced
            } else if (spanModel ||
                       ((intronModel || mixedModel) &&
                        !exonModel)) { // all models are only-spans. Not sure if this ever happens,
                                       // should be simplified. Intron or mixed, but not exon
                geneC[geneI][1]++; // unspliced
            } else { // all other combinations are mixed models, i.e. both only-exons and
                     // only-introns, only-exons and mixed
                geneC[geneI][2]++; // ambiguous
            };

            nUMIperCB[iCB]++;
        };

        countCellGeneUMIindex[iCB + 1] = countCellGeneUMIindex[iCB];

        if (nUMIperCB[iCB] == 0) // no UMIs counted for this CB
            continue;

        nGenePerCB[iCB] += geneC.size();
        readFeatSum->stats.V[readFeatSum->stats.yesUMIs] += nUMIperCB[iCB];
        ++readFeatSum->stats.V[readFeatSum->stats.yesCellBarcodes];

        if (countCellGeneUMI.size() <
            countCellGeneUMIindex[iCB + 1] +
                geneC.size() * countMatStride) // allocated vector too small
            countCellGeneUMI.resize(countCellGeneUMI.size() * 2);

        for (auto& gg : geneC) {
            countCellGeneUMI[countCellGeneUMIindex[iCB + 1] + 0] = gg.first;
            for (std::uint32_t ii = 0; ii < 3; ii++)
                countCellGeneUMI[countCellGeneUMIindex[iCB + 1] + 1 + ii] = gg.second[ii];
            countCellGeneUMIindex[iCB + 1] += countMatStride;
        };
    };

    // for now only unique mappers are considered for velocyto
    nReadPerCBtotal = nReadPerCB;
    nReadPerCBunique = nReadPerCB;

    time(&rawTime);
    context_.inOut->logMain << timeMonthDayTime(rawTime)
                            << " ... Velocyto counting: finished collapsing UMIs" << std::endl;
    context_.inOut->logMain << "RAM for solo feature " << FeatureTypes::Names[featureType] << "\n"
                            << linuxProcMemory() << std::flush;
};

} // namespace nustar::star2711b::solo::detail
