#include "../barcode/sequence.hpp"
#include "feature_counter.hpp"
#include "record_types.hpp"
#include "../support/algorithms.hpp"
#include "../support/streams.hpp"
#include "../support/time.hpp"
#include <bitset>
#include <math.h> /* log */
#include <numeric>
#include <unordered_map>

namespace nustar::star2711b::solo::detail {

void FeatureCounter::quantTranscript() { // velocyto counting gets info from Gene counting
    time_t rawTime;

    if (options.clusterCBfile == "-")
        return; // no transcript quantification w/o cluster file

    std::unordered_map<std::uint32_t, std::uint32_t>
        clusterCBind; // for each CB - cluster index of detected CB
    std::set<std::uint32_t>
        clusterInd; // cluster index for each cluster - the integer from clusterCBfile
    {               // load cluster information
        std::ifstream& clusterStream =
            ifstrOpen(options.clusterCBfile, sourceLocation(),
                      "SOLUTION: check the path and permissions of the cluster CB file: " +
                          options.clusterCBfile,
                      context_);
        std::string seq1;
        while (clusterStream >> seq1) {
            std::uint32_t icl1;
            clusterStream >> icl1;
            UnsignedWord cb1;
            if (convertNuclStrToInt64(seq1, cb1)) { // convert to 2-bit format
                auto cb1it = std::equal_range(options.cbWL.begin(), options.cbWL.end(),
                                              cb1); // find iterator in WL matching cb1
                std::uint32_t cb1ind =
                    (std::uint32_t)(cb1it.first -
                                    options.cbWL.begin()); // substract WL.begin iterator to find
                                                           // index in WL
                if (cb1ind < options.cbWL.size()) { // otherwise cb1 is not in WL
                    clusterCBind[cb1ind] = icl1;    // map: key=CB, value=cluster index
                    clusterInd.emplace(icl1);       // ordered set of cluster indexes
                } else {
                    context_.inOut->logMain
                        << "WARNING: cluster CB sequence not present in whitelist and is ignored: "
                        << seq1 << std::endl;
                };
            } else {
                context_.inOut->logMain
                    << "WARNING: cluster CB sequence contains non-ACGT base and is ignored: "
                    << seq1 << std::endl;
            };
        };
    };

    auto& trDistCount =
        readFeatSum->transcriptDistCount; // transcriptDistCount is accumulated while mapping, from
                                          // reads that map uniquely
    std::vector<double> trDistFun(trDistCount.size(), 0.0);
    std::vector<double> trDistFunTrFactor(annotation_.nTr, 0.0);

    { // process distance distribution function
        // running average
        std::int32_t runAverN = 50, runAverStart = 0;
        //    trDistFun[0] += trDistCount[ii];
        for (std::int32_t ii = 0; ii < runAverStart; ii++)
            trDistFun[ii] = (double)trDistCount[ii];
        for (std::int32_t ii = runAverStart; ii < (std::int32_t)trDistCount.size() - runAverN - 1;
             ii++) {
            trDistFun[ii] =
                (double)std::accumulate(trDistCount.begin() + std::max(runAverStart, ii - runAverN),
                                        trDistCount.begin() + ii + runAverN + 1, 0) /
                std::min(2 * runAverN + 1, ii - runAverStart + runAverN);
        };

        // cut when becomes non-monotonic
        std::uint32_t imax = 1000;
        while (trDistFun[imax + 1] > trDistFun[imax])
            imax++; // find maximum going forward from imax=1000
        context_.inOut->logMain << "SoloQuant: distance distribution past maximum = " << imax
                                << std::endl;

        while (trDistFun[imax + 1] < trDistFun[imax])
            imax++; // find first minimum after the maximum found above
        context_.inOut->logMain << "SoloQuant: distance distribution cutoff = " << imax
                                << std::endl;

        trDistFun.resize(imax);

        // normalize
        double norm1 = std::accumulate(trDistFun.begin(), trDistFun.end(), 0.0);

        std::ofstream* streamTrDistFun = &ofstrOpen(
            outputPrefix + "transcriptEndDistanceDistribution.txt", sourceLocation(), context_);
        for (auto& ff : trDistFun) {
            ff = ff / norm1;
            *streamTrDistFun << ff << '\n';
        };
        streamTrDistFun->close();

        // normalization factors for all transcripts
        std::vector<double> trDistFunCum(trDistFun);
        std::partial_sum(trDistFun.begin(), trDistFun.end(), trDistFunCum.begin());
        for (std::uint32_t ii = 0; ii < annotation_.nTr; ii++)
            if (annotation_.trLen[ii] < trDistFunCum.size())
                trDistFunTrFactor[ii] = -std::log(trDistFunCum[annotation_.trLen[ii] - 1]);

        for (auto& ff : trDistFun)
            ff = std::log(ff); // now trDistFun is log of dist function
    };

    typedef struct {
        std::uint32_t tr;
        double d;
    } transcriptDistProbStruct;
    // key=cluster index, key=CB/umi, value=vector of <trID, sum{trDistFun[distance3p]}> for each
    // read with the same umi
    std::map<std::uint32_t, std::unordered_map<UnsignedWord, std::vector<transcriptDistProbStruct>>>
        mapTrDist;
    for (int iThread = 0; iThread < context_.runThreadN;
         iThread++) { // TODO: this can be parallelized
        std::fstream* streamReads = readFeatAll[iThread]->streamReads;
        streamReads->flush();
        streamReads->seekg(0, std::ios::beg);

        std::string line1;
        while (std::getline(*streamReads, line1)) { // until the end of file
            std::stringstream lineStream(line1);
            std::uint32_t cb, cbCl, nTr;
            UnsignedWord umi;
            lineStream >> cb >> umi >> nTr;

            if (clusterCBind.count(cb) == 0) // this cb is not in the clusters
                continue;

            umi += (((UnsignedWord)cb) << 32); // to make UMI from different CB distinguishable

            cbCl = clusterCBind[cb];

            std::vector<transcriptDistProbStruct> tD;
            tD.reserve(nTr);
            for (std::uint32_t ii = 0; ii < nTr; ii++) {
                std::uint32_t tr1, d1;
                lineStream >> tr1 >> d1;
                if (d1 >= trDistFun.size())
                    continue; // do not record such outlier

                tD.push_back({tr1, trDistFun[d1] + trDistFunTrFactor[tr1]});
            };

            if (tD.size() == 0)
                continue;

            std::sort(tD.begin(), tD.end(),
                      [](const transcriptDistProbStruct& t1, const transcriptDistProbStruct& t2) {
                          return (t1.tr < t2.tr);
                      });

            if (mapTrDist[cbCl].count(umi) == 0) { // 1st entry for this umi
                mapTrDist[cbCl][umi] = tD;
                continue;
            };

            std::uint32_t inew = 0;
            std::vector<transcriptDistProbStruct> tD1;
            tD1.reserve(mapTrDist[cbCl][umi].size());

            for (std::uint32_t iold = 0; iold < mapTrDist[cbCl][umi].size();
                 iold++) { // intersection of old with new
                while (inew < tD.size() &&
                       mapTrDist[cbCl][umi][iold].tr > tD[inew].tr) // move through the sorted lists
                    ++inew; // inew advances if old_trID>new_trID

                if (inew == tD.size()) // end of tD reached
                    break;

                if (mapTrDist[cbCl][umi][iold].tr == tD[inew].tr) { // found match
                                                                    // old_trID==new_trID
                    tD1.push_back({tD[inew].tr,
                                   mapTrDist[cbCl][umi][iold].d +
                                       tD[inew].d}); // add new log probability to the existing one
                };
            };
            mapTrDist[cbCl][umi] = tD1; // replace with intersection
        };
    };

    time(&rawTime);
    context_.inOut->logMain << timeMonthDayTime(rawTime)
                            << " ... Transcript3p counting: finished input" << std::endl;

    std::map<std::uint32_t, std::vector<double>>
        clusterExpression; // per cluster, relative abundance
    // TODO parallelize this loop
    for (auto& mapTrDist1 : mapTrDist) { // loop over cell clusters
        auto& clTrDist = mapTrDist1.second;

        std::vector<double> trUnique(annotation_.nTr, 0),
            trInitial(annotation_.nTr, 0); // counts of unique read for each transcripts (i.e. reads
                                           // that map uniquely only to this transcript)
        UnsignedWord nUMItot = 0, nUMI0 = 0, nUMI1 = 0;

        { // pre-process mapTrDist:
            auto clTrDist1 = clTrDist.begin();
            while (clTrDist1 != clTrDist.end()) { // loop over UMIs
                auto& trDist = clTrDist1->second;

                if (trDist.size() == 0) {
                    clTrDist1 = clTrDist.erase(clTrDist1);
                    nUMI0++; // these are cases where the intersection of transcript from different
                             // reads with the same UMI is *empty*
                    continue;
                } else if (trDist.size() == 1) {
                    trUnique[trDist[0].tr]++;
                    trInitial[trDist[0].tr] += 1.0;
                    clTrDist1 = clTrDist.erase(clTrDist1);
                    nUMI1++; // only one transcript in the intersection, i.e. unique mappers
                    nUMItot++;
                    continue;
                };

                double max1 = std::max_element(
                                  trDist.begin(), trDist.end(),
                                  [](const transcriptDistProbStruct& t1,
                                     const transcriptDistProbStruct& t2) { return (t1.d < t2.d); })
                                  ->d;    // finds maximum d
                for (auto& tt : trDist) { // loop over transcripts in the UMI
                    trInitial[tt.tr] +=
                        1.0 / trDist.size(); // for initialization, split each count between
                                             // multimappers evenly and add it to unique mappers
                    tt.d = std::exp(tt.d - max1); // from sum(log(prob)) to product(prob), dividing
                                                  // by the max element to avoid underflow
                }; // this constant factor does not matter, since only the ratios of tt.d matter

                nUMItot++;
                clTrDist1++;
            };
            context_.inOut->logMain << timeMonthDayTime(rawTime)
                                    << " ... Transcript3p counting: cluster " << mapTrDist1.first
                                    << " nUMItot=" << nUMItot << " nUMI0=" << nUMI0
                                    << " nUMI1=" << nUMI1 << std::endl;
        };

        //         double trExpressed=0;
        //             trExpressed += tt;
        //
        //             tt *= nUMItot/trExpressed;

        std::vector<double> thOldNew[2];
        thOldNew[0] = trInitial;
        thOldNew[1].resize(annotation_.nTr, 0);

        auto* thOldP = &thOldNew[0];
        auto* thNewP = &thOldNew[1];

        std::vector<bool> trConverged(
            trInitial.size(),
            false); // if true, this transcript is converged - then it's not updated in EM

        for (std::uint32_t iteration = 0; iteration < 10000; iteration++) { // main EM loop
            auto& thNew = *thNewP;
            auto& thOld = *thOldP;

            // always start with trUnique counts
            std::copy(trUnique.begin(), trUnique.end(), thNew.begin());

            // calculate thNew
            for (auto& clTrDist1 : clTrDist) { // loop over all UMIs
                auto& trDist = clTrDist1.second;

                double denom1 = 0.0;
                for (auto& td : trDist) { // loop over transcripts in one UMI: caclulate denomiator
                    denom1 += td.d * thOld[td.tr];
                };

                for (auto& td : trDist) { // loop over transcripts in one UMI: update thNew
                    if (!trConverged[td.tr])
                        thNew[td.tr] += td.d * thOld[td.tr] / denom1; // adding to unique counts
                };
            };

            // convergence check
            double diffThresholdMax = 1e-5;
            double diffThresholdOne = diffThresholdMax * 0.1;
            double exprThreshold = 1e-8 * nUMItot;
            double diffMax = 0, diffSum = 0, aboveThrN = 0, aboveThrExprSum = 0, aboveThrOneN = 0;
            for (std::uint32_t itr = 0; itr < thNew.size(); itr++) {
                if (trConverged[itr] || thOld[itr] == 0)
                    continue;

                double diff1 = std::abs(thNew[itr] - thOld[itr]) / thOld[itr];
                diffSum += diff1;
                diffMax = std::max(diffMax, diff1);
                if (diff1 > diffThresholdMax) {
                    aboveThrN++;
                    aboveThrExprSum += thNew[itr];
                };

                if (thNew[itr] < exprThreshold) {
                    trConverged[itr] = true;
                    // trUnique[itr]=thNew[itr];//0 here could make more sense, but convergence
                    // suffers
                    trUnique[itr] = 0;
                };
                if (diff1 < diffThresholdOne) {
                    trConverged[itr] = true;
                    trUnique[itr] = thNew[itr]; // this is the final value for this transcript
                } else {
                    aboveThrOneN++;
                };
            };

            std::cout << iteration << " " << diffMax << " " << diffSum << " " << aboveThrN << " "
                      << aboveThrExprSum << " " << aboveThrOneN << std::endl;

            if (diffMax < diffThresholdMax)
                break;
            std::swap(thNewP, thOldP); // swap for the next iteration
        };

        auto& thOut = *thNewP;
        { // renormalize into TPMs
            double norm1 = 0;
            for (std::uint32_t itr = 0; itr < thOut.size(); itr++) {
                thOut[itr] *= std::exp(trDistFunTrFactor[itr]);
                norm1 += thOut[itr];
            };
            norm1 = nUMItot / norm1;
            for (std::uint32_t itr = 0; itr < thOut.size(); itr++) {
                thOut[itr] *= norm1;
            };
        };

        clusterExpression[mapTrDist1.first] = thOut;

        time(&rawTime);
        context_.inOut->logMain << timeMonthDayTime(rawTime)
                                << " ... Transcript3p counting: finished cluster"
                                << mapTrDist1.first << std::endl;
    };

    { // output counting matrix
        std::string matrixFileName = outputPrefix + options.outFileNames[3];
        std::ofstream& countMatrixStream = ofstrOpen(matrixFileName, sourceLocation(), context_);
        countMatrixStream << "%%MatrixMarket matrix coordinate real general\n%\n";

        std::uint32_t nCellGeneEntries = 0;
        for (auto& ctpm : clusterExpression)
            for (auto& tt : ctpm.second)
                if (tt > 0)
                    nCellGeneEntries++;

        countMatrixStream << annotation_.nTr << ' '
                          << *std::max_element(clusterInd.begin(), clusterInd.end()) << ' '
                          << nCellGeneEntries << '\n';

        for (auto& ctpm : clusterExpression)
            for (auto itt = ctpm.second.begin(); itt != ctpm.second.end(); itt++)
                if (*itt > 0)
                    countMatrixStream << itt - ctpm.second.begin() + 1 << ' ' << ctpm.first << ' '
                                      << *itt << '\n';

        countMatrixStream.flush();
    };

    std::ofstream& outStr = ofstrOpen(outputPrefix + "/features.tsv", sourceLocation(), context_);
    for (std::uint32_t ii = 0; ii < annotation_.nTr; ii++)
        outStr << annotation_.trID[ii] << "\t" << annotation_.trLen[ii] << "\t"
               << annotation_.geName[annotation_.trGene[ii]] << '\n';
    outStr.close();

    time(&rawTime);
    context_.inOut->logMain << timeMonthDayTime(rawTime)
                            << " ... Transcript3p counting: finished transcript quantification"
                            << std::endl;
};

} // namespace nustar::star2711b::solo::detail
