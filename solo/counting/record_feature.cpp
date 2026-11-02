#include "../annotation/read_annotations.hpp"
#include "../barcode/read_barcode.hpp"
#include "../barcode/sequence.hpp"
#include "feature_recorder.hpp"
#include "../support/algorithms.hpp"

namespace nustar::star2711b::solo::detail {

class ReadSoloFeatures {
public:
    std::uint32_t gene;
    std::vector<std::uint32_t> geneMult;
    std::vector<std::array<UnsignedWord, 2>> sj;
    bool sjAnnot;
    std::uint32_t indAnnotTr; // index of the annotated transcript
    Alignment** alignOut;
};

std::uint32_t outputReadCB(std::fstream* streamOut, const UnsignedWord iRead,
                           const std::int32_t featureType, ReadBarcode& soloBar,
                           const ReadSoloFeatures& reFe, const ReadAnnotations& readAnnot,
                           const ReadFlags& readFlag);

void FeatureRecorder::record(ReadBarcode& soloBar, UnsignedWord nTr, Alignment** alignOut,
                             UnsignedWord iRead, ReadAnnotations& readAnnot) {
    if (options.type == 0)
        return;

    if (options.readStatsYes[featureType]) { // readFlag

        if (nTr == 1) {
            readFlag.setBit(readFlag.genomeU);
        } else if (nTr > 1) {
            readFlag.setBit(readFlag.genomeM);
        };

        for (UnsignedWord itr = 0; itr < nTr; itr++) {
            if (context_.pGe.chrSet.mito.count(alignOut[itr]->Chr) == 1) {
                readFlag.setBit(readFlag.mito);
            };
        };

        switch (readAnnot.annotFeatures[featureType].ovType) {
        case ReadAnnotFeature::overlapTypes::exonic:
        case ReadAnnotFeature::overlapTypes::exonic50p:
            readFlag.setBit(readFlag.exonic);
            break;
        case ReadAnnotFeature::overlapTypes::intronic:
            readFlag.setBit(readFlag.intronic);
            break;
        case ReadAnnotFeature::overlapTypes::exonicAS:
        case ReadAnnotFeature::overlapTypes::exonic50pAS:
            readFlag.setBit(readFlag.exonicAS);
            break;
        case ReadAnnotFeature::overlapTypes::intronicAS:
            readFlag.setBit(readFlag.intronicAS);
        };

        if (soloBar.cbMatch < 0 && options.cbWLyes) { // no CB match in the WL
            if (readAnnot.annotFeatures[featureType].fSet.size() == 1) {
                readFlag.setBit(readFlag.featureU);
            } else if (readAnnot.annotFeatures[featureType].fSet.size() > 1) {
                readFlag.setBit(readFlag.featureM);
            };
            readFlag.setBit(readFlag.cbMatch); // this will counts reads with no CB match
            readFlag.countsAddNoCB();
        };
    };

    if (soloBar.cbMatch < 0)
        return;

    ReadSoloFeatures reFe;
    reFe.alignOut = alignOut;
    reFe.indAnnotTr = 0;

    std::uint32_t nFeat = 0; // number of features in this read (could be >1 for SJs)
    if (nTr == 0) {          // unmapped
        stats.V[stats.noUnmapped]++;

    } else {
        switch (featureType) {
        case FeatureTypes::Gene:
        case FeatureTypes::GeneFull:
        case FeatureTypes::GeneFull_Ex50pAS:
        case FeatureTypes::GeneFull_ExonOverIntron: {
            auto* readGe = &readAnnot.annotFeatures[featureType].fSet;

            if (soloBar.options.type == soloBar.options.SoloTypes::SmartSeq) {
                for (std::int32_t itr = nTr - 1; itr >= 0; itr--) {
                    if (readAnnot.annotFeatures[featureType].fAlign[itr].size() > 0) {
                        reFe.indAnnotTr = itr;
                        break; // indAnnotTr is the first genic align, and is used for Smart-seq
                               // deduplication, below in outputReadCB()
                               // TODO this is ambiguous for reads that map to the same gene
                               // multiple times. *Last* align is chosen to match the previous
                               // version.
                    };
                };
            };

            if (readGe->size() == 0) {        // check genes
                stats.V[stats.noNoFeature]++; // no gene
            } else if (readGe->size() > 1) {
                stats.V[stats.MultiFeature]++; // multi-gene reads
                readFlag.setBit(readFlag.featureM);
                if (nTr > 1)
                    stats.V[stats.subMultiFeatureMultiGenomic]++; // multigene caused by multimapper

                if (options.multiMap.yes.multi) { // output multimappers
                    reFe.geneMult.resize(readGe->size());
                    std::uint32_t ii = 0;
                    for (auto& g : *readGe) { // set high bit to mark multimappers
                        reFe.geneMult[ii] = g | geneMultMark;
                        ++ii;
                    };

                    nFeat = outputReadCB(streamReads, iRead, featureType, soloBar, reFe, readAnnot,
                                         readFlag);
                };
            } else { // unique-gene reads
                reFe.gene = *readGe->begin();
                readFlag.setBit(readFlag.featureU);
                nFeat = outputReadCB(streamReads, (readIndexYes ? iRead : (UnsignedWord)-1),
                                     featureType, soloBar, reFe, readAnnot, readFlag);
            };

            // debug
            //{
            //         cout << iRead<<" "<<readFlag.checkBit(readFlag.featureM)<<"
            //         "<<readFlag.checkBit(readFlag.featureU)<<"
            //         "<<readFlag.checkBit(readFlag.exonic)<<endl;
            // };
        }; break;

        case FeatureTypes::SJ:
            if (nTr > 1) { // reject all multimapping reads
                stats.V[stats.subMultiFeatureMultiGenomic]++;
                stats.V[stats.MultiFeature]++;
                //} else if (readAnnot.geneConcordant.size()>1){//for SJs, still check genes, no
                //feature if multi-gene
                //    stats.V[stats.MultiFeature]++;
            } else { // one gene or no gene
                alignOut[0]->extractSpliceJunctions(reFe.sj, reFe.sjAnnot);
                if (reFe.sj.empty()) {
                    stats.V[stats.noNoFeature]++;
                } else { // good junction
                    readFlag.setBit(readFlag.featureU);
                    nFeat = outputReadCB(streamReads, (readIndexYes ? iRead : (UnsignedWord)-1),
                                         featureType, soloBar, reFe, readAnnot, readFlag);
                };
            };
            break;

        case FeatureTypes::Transcript3p:
            if (readAnnot.transcriptConcordant.size() == 0 ||
                soloBar.cbMatch > 1) { // do not record ambiguous CB
                stats.V[stats.noNoFeature]++;
            } else {
                nFeat = outputReadCB(streamReads, iRead, featureType, soloBar, reFe, readAnnot,
                                     readFlag);
            };
            if (readAnnot.transcriptConcordant.size() == 1 &&
                readAnnot.transcriptConcordant[0][1] < transcriptDistCount.size()) {
                // read maps to one transcript - use for distTTS distribution function
                transcriptDistCount[readAnnot.transcriptConcordant[0][1]]++;
            };

            break;

        case FeatureTypes::Velocyto:
            // different record: iRead, nTr, tr1, type1, tr2, type2 ...
            if (readAnnot.trVelocytoType.size() > 0) { // otherwise, no gene

                std::sort(readAnnot.trVelocytoType.begin(), readAnnot.trVelocytoType.end(),
                          [](const TranscriptType& t1, const TranscriptType& t2) {
                              return t1.tr < t2.tr;
                          });

                *streamReads << iRead << ' ' << readAnnot.trVelocytoType.size();
                for (auto& tt : readAnnot.trVelocytoType)
                    *streamReads << ' ' << tt.tr << ' ' << (std::uint32_t)tt.type;
                *streamReads << '\n';
                nFeat = 1;
            } else {
                stats.V[stats.noNoFeature]++;
            };
            break; // no need to go with downstream processing

        }; // switch (featureType)
    }; // if (nTr==0)

    if (nFeat == 0 &&
        (readInfoYes | options.readStatsYes[featureType])) { // no feature, but readInfo requested
        outputReadCB(streamReads, iRead, (std::uint32_t)-1, soloBar, reFe, readAnnot, readFlag);
    };

    if (nFeat == 0)
        return; // no need to record the number of reads per CB

    if (options.cbWLyes) { // WL
        for (auto& cbi : soloBar.cbMatchInd)
            cbReadCount[cbi] += nFeat;
    } else { // no WL
        cbReadCountMap[soloBar.cbMatchInd[0]] += nFeat;
    };

    /*//if we wanted to record all barcodes, even those without features
        if (!pSolo.cbWLyes) {//no WL
        cbReadCountMap[soloBar.cbMatchInd[0]] += nFeat;
    } else if (nFeat>0){//WL && nFeat>0
        for (auto &cbi : soloBar.cbMatchInd)
            cbReadCount[cbi] += nFeat;
    };
    */

    return;
};

std::uint32_t outputReadCB(std::fstream* streamOut, const UnsignedWord iRead,
                           const std::int32_t featureType, ReadBarcode& soloBar,
                           const ReadSoloFeatures& reFe, const ReadAnnotations& readAnnot,
                           const ReadFlags& readFlag) {
    /*format of the temp output file
     * UMI [iRead] type feature* cbMatchString
     *             0=exact match, 1=one non-exact match, 2=multipe non-exact matches
     *                   gene or sj[0] sj[1]
     *                         CB or nCB {CB Qual, ...}
     */

    if (soloBar.options.type == soloBar.options.SoloTypes::SmartSeq &&
        featureType != -1) { // need to calculate "UMI" from align start/end
        soloBar.umiB = reFe.alignOut[reFe.indAnnotTr]->chrStartLengthExtended();
    };

    UnsignedWord nout = 1;

    switch (featureType) {
    case -1:
        // no feature, output for readInfo
        *streamOut << soloBar.umiB << ' ' << iRead << ' ' << readFlag.flag << ' ' << -1 << ' '
                   << soloBar.cbMatch << ' ' << soloBar.cbMatchString << '\n';
        break;

    case FeatureTypes::Gene:
    case FeatureTypes::GeneFull:
    case FeatureTypes::GeneFull_Ex50pAS:
    case FeatureTypes::GeneFull_ExonOverIntron:
        if (reFe.geneMult.size() == 0) {
            // just gene id
            *streamOut << soloBar.umiB << ' '; // UMI
            if (iRead != (UnsignedWord)-1)
                *streamOut << iRead << ' ' << readFlag.flag << ' '; // iRead
            *streamOut << reFe.gene << ' ' << soloBar.cbMatch << ' ' << soloBar.cbMatchString
                       << '\n';
            ;
        } else {
            for (auto& g : reFe.geneMult) {
                *streamOut << soloBar.umiB << ' '; // UMI
                *streamOut << iRead << ' ' << readFlag.flag
                           << ' '; // iRead is always output for multiGene
                *streamOut << g << ' ' << soloBar.cbMatch << ' ' << soloBar.cbMatchString << '\n';
                ;
            };
            nout = reFe.geneMult.size();
        };
        break;

    case FeatureTypes::SJ:
        // sj - two numbers, multiple sjs per read
        for (auto& sj : reFe.sj) {
            *streamOut << soloBar.umiB << ' '; // UMI
            if (iRead != (UnsignedWord)-1)
                *streamOut << iRead << ' ' << readFlag.flag << ' '; // iRead
            *streamOut << sj[0] << ' ' << sj[1] << ' ' << soloBar.cbMatch << ' '
                       << soloBar.cbMatchString << '\n'
                       << std::flush;
        };
        nout = reFe.sj.size();
        break;

    case FeatureTypes::Transcript3p:
        // transcript,distToTTS structure
        *streamOut << soloBar.cbMatchString << ' ';
        *streamOut << soloBar.umiB << ' ';
        *streamOut << readAnnot.transcriptConcordant.size();
        for (auto& tt : readAnnot.transcriptConcordant) {
            *streamOut << ' ' << tt[0] << ' ' << tt[1];
        };
        if (iRead != (UnsignedWord)-1)
            *streamOut << ' ' << iRead; // iRead
        *streamOut << '\n';
        nout = 1;

        break;
    }; // switch (featureType)

    return nout;
};

} // namespace nustar::star2711b::solo::detail
