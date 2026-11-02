#include "../barcode/sequence.hpp"
#include "feature_counter.hpp"
#include "../support/streams.hpp"
#include "../support/time.hpp"

namespace nustar::star2711b::solo::detail {

void FeatureCounter::sumThreads() {
    // stats
    nReadsInput = context_.countingStats.readN + 1;    // reserve 1 extra
    for (int ii = 0; ii < context_.runThreadN; ii++) { // point to
        readFeatAll[ii] = workers[ii]->read->readFeat[options.featureInd[featureType]].get();
        readFeatAll[ii]->streamReads->flush();
        readFeatSum->addCounts(*readFeatAll[ii]);
    };

    // if WL was not defined
    if (!options.cbWLyes) { // now we can define WL and counts ??? we do not need to do it for every
                            // feature???
        options.cbWLsize = readFeatSum->cbReadCountMap.size();
        options.cbWL.resize(options.cbWLsize);
        options.cbWLstr.resize(options.cbWLsize);
        readFeatSum->cbReadCount.resize(options.cbWLsize);
        barcode_sum_->cbReadCountExact.resize(options.cbWLsize);

        if (options.CBtype.type == 1) { // sequence cb
            UnsignedWord icb = 0;
            for (auto& cb : readFeatSum->cbReadCountMap) {
                options.cbWL[icb] = cb.first;
                options.cbWLstr[icb] = convertNuclInt64toString(options.cbWL[icb], options.cbL);
                readFeatSum->cbReadCount[icb] = cb.second;
                barcode_sum_->cbReadCountExact[icb] = cb.second;
                icb++;
            };
        } else if (options.CBtype.type == 2) { // string cb
            std::vector<std::unordered_map<std::string, std::uint32_t>::iterator> cbiter(
                options.CBtype.strMap.size());
            for (auto cbi = options.CBtype.strMap.begin(); cbi != options.CBtype.strMap.end();
                 cbi++)
                cbiter[cbi->second] = cbi;

            UnsignedWord icb = 0;
            for (auto& cb : readFeatSum->cbReadCountMap) {
                options.cbWL[icb] = cb.first;
                options.cbWLstr[icb] = cbiter[cb.first]->first;
                readFeatSum->cbReadCount[icb] = cb.second;
                barcode_sum_->cbReadCountExact[icb] = cb.second;
                icb++;
            };
        };

        // pseudocounts
        if (options.CBmatchWL.mm1_multi_pc) {
            for (std::uint32_t ii = 0; ii < options.cbWLsize; ii++) {
                barcode_sum_->cbReadCountExact[ii]++; // add one to exact counts
            };
        };
    };

    // if restarting from _STARtmp/solo* file
    if (context_.runRestart.type == 1) { // this could happen if the run is restarted. Would be
                                         // better to save/load cbReadCount, or recalculate it from
        for (int ii = 0; ii < context_.runThreadN; ii++) {
            readFeatAll[ii]
                ->streamReads->clear(); // just in case EOF was reached in previous reading
            readFeatAll[ii]->streamReads->seekg(0, std::ios::beg);
            std::string line1;
            while (std::getline(*readFeatAll[ii]->streamReads, line1)) {
                std::istringstream line1stream(line1);
                UnsignedWord cb1;
                line1stream >> cb1 >> cb1 >> cb1;
                if (featureType == FeatureTypes::SJ)
                    line1stream >> cb1;
                line1stream >> cb1;
                //    continue;//this should not happen!
                readFeatSum->cbReadCount[cb1]++;
            };
        };
    };

    // detected CBs
    nCB = 0;
    nReadsMapped = 0;
    for (std::uint32_t ii = 0; ii < options.cbWLsize; ii++) {
        if (readFeatSum->cbReadCount[ii] > 0) {
            nCB++;
            nReadsMapped += readFeatSum->cbReadCount[ii];
        };
    };

    indCBwl.resize(options.cbWLsize, (std::uint32_t)-1);
    indCB.resize(nCB);
    nCB = 0; // will count it again below
    for (std::uint32_t ii = 0; ii < options.cbWLsize; ii++) {
        if (readFeatSum->cbReadCount[ii] > 0) {
            indCB[nCB] = ii;
            indCBwl[ii] = nCB;
            ++nCB;
        };
    };
};

} // namespace nustar::star2711b::solo::detail
