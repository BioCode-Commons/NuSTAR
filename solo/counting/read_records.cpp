#include "feature_counter.hpp"
#include "feature_recorder.hpp"
#include "input_record.hpp"
#include "record_types.hpp"
#include "../support/algorithms.hpp"
#include <cmath>

namespace nustar::star2711b::solo::detail {

void FeatureRecorder::inputRecords(std::uint32_t** cbP, std::uint32_t cbPstride,
                                   std::vector<std::uint32_t>& cbReadCountTotal,
                                   std::vector<ReadCorrection>& readInfo, ReadFlags& readFlagCounts,
                                   std::vector<std::uint32_t>& nReadPerCBunique1,
                                   std::vector<std::uint32_t>& nReadPerCBmulti1) {
    streamReads->flush();
    streamReads->seekg(0, std::ios::beg);
    std::uint32_t feature;
    UnsignedWord umi, iread, prevIread = (UnsignedWord)-1;
    std::int32_t cbmatch;
    SignedWord cb;
    std::vector<std::uint32_t> trIdDist;

    UnsignedWord nReadsIn = 0;

    while (soloInputFeatureUMI(streamReads, featureType, readIndexYes, context_.sjAll, iread,
                               cbmatch, feature, umi, trIdDist, readFlagCounts)) {
        if (feature == (std::uint32_t)(-1) &&
            !readIndexYes) { // no feature => no record, this can happen for SJs
            streamReads->ignore((std::uint32_t)-1, '\n');
            // stats.V[stats.noNoFeature]++; //need separate category for this
            continue;
        };

        bool readIsCounted = false;
        bool featGood = (feature != (std::uint32_t)(-1));
        bool noMMtoWLwithoutExact = false;
        bool noTooManyWLmatches = false;

        if (cbmatch <= 1) { // single match
            *streamReads >> cb;

            if (options.CBmatchWL.oneExact && cbmatch == 1 &&
                cbReadCountTotal[cb] == 0) { // single 1MM match, no exact matches to this CB
                noMMtoWLwithoutExact = true;

            } else {

                if (!options.cbWLyes) { // if no-WL, the full cbInteger was recorded - now has to be
                                        // placed in order
                    cb = binarySearchExact<BarcodeCode>(cb, options.cbWL.data(), options.cbWLsize);
                    if (cb + 1 == 0)
                        continue; // this cb was not in the tentative WL
                };

                // record feature
                if (featGood) { // good feature, will be counted
                    readIsCounted = true;

                    cbP[cb][0] = feature;
                    cbP[cb][1] = umi;

                    if (readIndexYes) {
                        cbP[cb][2] = iread;
                    };

                    cbP[cb] += cbPstride;

                } else if (readInfoYes) { // no feature - record readInfo
                    readInfo[iread].cb = cb;
                    readInfo[iread].umi = umi;
                };
            };

        } else { // multiple matches
#ifdef MATCH_CellRanger
            double ptot = 0.0, pmax = 0.0, pin;
#else
            float ptot = 0.0, pmax = 0.0, pin;
#endif

            for (std::uint32_t ii = 0; ii < (std::uint32_t)cbmatch; ii++) {
                std::uint32_t cbin;
                char qin;
                *streamReads >> cbin >> qin;
                if (cbReadCountTotal[cbin] > 0) { // otherwise this cbin does not work
                    qin -= options.QSbase;
                    qin = qin < options.QSmax ? qin : options.QSmax;
                    pin = cbReadCountTotal[cbin] * std::pow(10.0, -qin / 10.0);
                    ptot += pin;
                    if (pin > pmax) {
                        cb = cbin;
                        pmax = pin;
                    };
                };
            };
            if (ptot > 0.0 && pmax >= options.cbMinP * ptot) {
                // record feature single-number feature
                if (featGood) {
                    readIsCounted = true;
                    cbP[cb][0] = feature;
                    cbP[cb][1] = umi;
                    if (readIndexYes) {
                        cbP[cb][2] = iread;
                    };
                    cbP[cb] += cbPstride;
                } else if (readInfoYes) { // no feature - record readInfo
                    readInfo[iread].cb = cb;
                    readInfo[iread].umi = umi;
                };
            } else {
                noTooManyWLmatches = true;
            };
        };

        if (!readIndexYes ||
            iread != prevIread) { // no readindex (then only unique-gene reads are recorded) OR only
                                  // for one align of each read, in case of multimappers
            prevIread = iread;
            if (featGood) {
                if (cbmatch == 0) {
                    stats.V[stats.yessubWLmatchExact]++;
                } else if (noMMtoWLwithoutExact) {
                    stats.V[stats.noMMtoWLwithoutExact]++;
                } else if (noTooManyWLmatches) {
                    stats.V[stats.noTooManyWLmatches]++;
                };
            };

            if (readIsCounted) {
                if (feature < geneMultMark) {
                    nReadPerCBunique1[cb]++;
                } else {
                    nReadPerCBmulti1[cb]++;
                };
            };

            if (options.readStatsYes[featureType]) { // has to be new iread to avoid muti-counting
                                                     // multi-gene reads

                // readIsCounted flag was defined above
                if (readIsCounted) {
                    if (readFlagCounts.checkBit(readFlagCounts.featureU))
                        readFlagCounts.setBit(readFlagCounts.countedU);
                    if (readFlagCounts.checkBit(readFlagCounts.featureM))
                        readFlagCounts.setBit(readFlagCounts.countedM);
                };

                readFlagCounts.setBit(readFlag.cbMatch); // set this flag for both CB and no-CB, but
                                                         // no-CB will be counted separately below
                // only record this read in readFlagCounts if its CB was defined
                if (cbmatch == 0) { // perfect CB match to WL
                    readFlagCounts.setBit(readFlagCounts.cbPerfect);
                    readFlagCounts.countsAdd(cb);
                } else if (cbmatch == 1 && !noMMtoWLwithoutExact) {
                    readFlagCounts.setBit(readFlagCounts.cbMMunique);
                    readFlagCounts.countsAdd(cb);
                } else if (cbmatch > 1 && !noTooManyWLmatches) {
                    readFlagCounts.setBit(readFlagCounts.cbMMmultiple);
                    readFlagCounts.countsAdd(cb);
                } else { // no CB match
                    readFlagCounts.countsAddNoCB();
                };

                // debug
                /*nReadsIn++;
                uint64 ntot=0;
                for (auto &ii: readFlagCounts.flagCounts)
                    ntot += ii.second[0];

                if (nReadsIn!=ntot)
                    cout << nReadsIn <<" "<< ntot << endl;
                */
                /*//if (featureType==SoloFeatureTypes::Gene &&
                readFlagCounts.flagCounts[cb][readFlagCounts.featureM]+readFlagCounts.flagCounts[cb][readFlagCounts.featureU]
                != readFlagCounts.flagCounts[cb][readFlagCounts.exonic])
                //    cout << cb <<' '<< iread << endl;
                //    cout << iread <<' '<< readFlagCounts.flagCounts[cb][readFlagCounts.featureM]
                << endl;
                */
            };
        };
    };
};

} // namespace nustar::star2711b::solo::detail
