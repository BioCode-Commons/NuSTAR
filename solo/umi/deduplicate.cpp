#include "../counting/feature_counter.hpp"
#include "../counting/record_types.hpp"
#include "../support/algorithms.hpp"
#include "../support/streams.hpp"
#include "../support/time.hpp"
#include <unordered_map>

namespace nustar::star2711b::solo::detail {

inline int funCompareSolo1(const void* a, const void* b); // defined below
inline int funCompare_uint32_1_2_0(const void* a, const void* b);

void FeatureCounter::deduplicateAllCells() {
    std::vector<std::uint32_t> umiArray(nReadPerCBmax *
                                        umiArrayStride); // temp array for collapsing UMI
    std::vector<std::uint32_t> gID(
        std::min(2 * featuresNumber, nReadPerCBmax) +
        1); // gene IDs, 2* is needed because each gene can have unique and multi-mappers
    std::vector<std::uint32_t> gReadS(
        std::min(2 * featuresNumber, nReadPerCBmax) +
        1); // start of gene reads TODO: allocate this array in the 2nd half of rGU

    for (std::uint32_t icb = 0; icb < nCB; icb++) { // main collapse cycle

        deduplicateCell(icb, umiArray, gID, gReadS);

        readFeatSum->stats.V[readFeatSum->stats.yesUMIs] += nUMIperCB[icb];
        if (nGenePerCB[icb] > 0) // nGenePerCB contains only unique
            ++readFeatSum->stats.V[readFeatSum->stats.yesCellBarcodes];

        readFeatSum->stats.V[readFeatSum->stats.yesWLmatch] += nReadPerCBtotal[icb];
        readFeatSum->stats.V[readFeatSum->stats.yessubWLmatch_UniqueFeature] +=
            nReadPerCBunique[icb];
    };
};

void FeatureCounter::deduplicateCell(std::uint32_t iCB, std::vector<std::uint32_t>& umiArray,
                                     std::vector<std::uint32_t>& gID,
                                     std::vector<std::uint32_t>& gReadS) {

    std::uint32_t* rGU = rCBp[iCB];
    std::uint32_t rN =
        nReadPerCB[iCB]; // with multimappers, this is the number of all aligns, not reads

    qsort(rGU, rN, rguStride * sizeof(std::uint32_t),
          funCompareNumbers<std::uint32_t>); // sort by gene index

    std::uint32_t gid1 = -1;                  // current gID
    std::uint32_t nGenes = 0, nGenesMult = 0; // number of genes
    for (std::uint32_t iR = 0; iR < rN * rguStride; iR += rguStride) {
        if (rGU[iR + rguG] != gid1) { // record gene boundary
            gReadS[nGenes] = iR;
            gid1 = rGU[iR + rguG];
            gID[nGenes] = gid1;

            ++nGenes;
            if (options.multiMap.yes.multi && (gid1 & geneMultMark))
                ++nGenesMult;
        };
    };
    gReadS[nGenes] = rguStride * rN; // so that gReadS[nGenes]-gReadS[nGenes-1] is the number of
                                     // reads for nGenes, see below in qsort
    nGenes -= nGenesMult; // unique only gene

    // unordered_map <uintUMI, unordered_set<uint32>> umiGeneMap;
    std::unordered_map<UmiCode, std::unordered_map<std::uint32_t, std::uint32_t>> umiGeneMapCount,
        umiGeneMapCount0;
    // UMI                 //Gene //Count
    if (options.umiFiltering.MultiGeneUMI) {
        for (std::uint32_t iR = 0; iR < gReadS[nGenes]; iR += rguStride) {
            umiGeneMapCount[rGU[iR + 1]][rGU[iR]]++;
        };

        for (auto& iu : umiGeneMapCount) { // loop over all UMIs
            if (iu.second.size() == 1)
                continue;
            std::uint32_t maxu = 0;
            for (const auto& ig : iu.second) { // loop over genes for a given UMI
                if (maxu < ig.second)
                    maxu = ig.second; // find gene with maximum count
            };
            if (maxu == 1)
                maxu = 2; // to kill UMIs with 1 read to one gene, 1 read to another gene
            for (auto& ig : iu.second) {
                if (maxu > ig.second)
                    ig.second = 0; // kills Gene with read count *strictly* < maximum count
            };
        };
    };

    if (options.umiFiltering.MultiGeneUMI_All) {
        for (std::uint32_t iR = 0; iR < gReadS[nGenes]; iR += rguStride) {
            umiGeneMapCount[rGU[iR + 1]][rGU[iR]]++;
        };

        for (auto& iu : umiGeneMapCount) { // loop over all UMIs
            if (iu.second.size() > 1) {
                for (auto& ig : iu.second)
                    ig.second = 0; // kill all genes for this UMI
            };
        };
    };

    std::vector<std::unordered_map<UmiCode, UmiCode>> umiCorrected(nGenes);

    if (countCellGeneUMI.size() < countCellGeneUMIindex[iCB] + nGenes * countMatStride)
        countCellGeneUMI.resize((countCellGeneUMI.size() + nGenes * countMatStride) *
                                2); // allocated vector too small

    nGenePerCB[iCB] = 0;
    nUMIperCB[iCB] = 0;
    countCellGeneUMIindex[iCB + 1] = countCellGeneUMIindex[iCB];
    for (std::uint32_t iG = 0; iG < nGenes; iG++) { // collapse UMIs for each gene
        std::uint32_t* rGU1 = rGU + gReadS[iG];

        std::uint32_t nR0 = (gReadS[iG + 1] - gReadS[iG]) / rguStride; // total number of reads
        if (nR0 == 0)
            continue; // no reads - this should not happen?

        qsort(rGU1, nR0, rguStride * sizeof(std::uint32_t),
              funCompareTypeShift<std::uint32_t, rguU>);

        // exact collapse
        std::uint32_t iR1 = -umiArrayStride; // number of distinct UMIs for this gene
        std::uint32_t u1 = -1;
        for (std::uint32_t iR = rguU; iR < gReadS[iG + 1] - gReadS[iG];
             iR += rguStride) { // count and collapse identical UMIs
            if (options.umiFiltering.MultiGeneUMI &&
                umiGeneMapCount[rGU1[iR]][gID[iG]] == 0) { // multigene UMI is not recorded
                if (options.umiDedup.typeMain !=
                    UmiDedupOptions::typeI::NoDedup) // for NoDedup, the UMI filtering is not done
                    rGU1[iR] = (UmiCode)-1; // mark multigene UMI, so that UB tag will be set to -
                continue;
            };

            if (rGU1[iR] != u1) {
                iR1 += umiArrayStride;
                u1 = rGU1[iR];
                umiArray[iR1] = u1;
                umiArray[iR1 + 1] = 0;
            };
            umiArray[iR1 + 1]++;
        };

        std::uint32_t nU0 =
            (iR1 + umiArrayStride) / umiArrayStride; // number of UMIs after simple exact collapse

        if (options.umiFiltering.MultiGeneUMI_CR) {
            if (nU0 == 0)
                continue; // nothing to count

            for (UnsignedWord iu = 0; iu < nU0 * umiArrayStride; iu += umiArrayStride) {
                umiGeneMapCount0[umiArray[iu + 0]][iG] +=
                    umiArray[iu + 1]; // this sums read counts over UMIs that were collapsed
            };

            correctUmisCellRanger(nU0, umiArray.data(), readInfo.size() > 0, false,
                                  umiCorrected[iG]);

            for (UnsignedWord iu = 0; iu < nU0 * umiArrayStride;
                 iu += umiArrayStride) { // just fill the umiGeneMapCount - will calculate UMI
                                         // counts later
                umiGeneMapCount[umiArray[iu + 2]][iG] +=
                    umiArray[iu + 1]; // this sums read counts over UMIs that were collapsed
            };

            continue; // done with MultiGeneUMI_CR, readInfo will be filled later
        };

        if (options.umiDedup.yes.NoDedup)
            countCellGeneUMI[countCellGeneUMIindex[iCB + 1] + options.umiDedup.countInd.NoDedup] =
                nR0;

        if (nU0 > 0) { // otherwise no need to count
            if (options.umiDedup.yes.Exact)
                countCellGeneUMI[countCellGeneUMIindex[iCB + 1] + options.umiDedup.countInd.Exact] =
                    nU0;

            if (options.umiDedup.yes.CR)
                countCellGeneUMI[countCellGeneUMIindex[iCB + 1] + options.umiDedup.countInd.CR] =
                    correctUmisCellRanger(nU0, umiArray.data(),
                                          readInfo.size() > 0 && options.umiDedup.typeMain ==
                                                                     UmiDedupOptions::typeI::CR,
                                          true, umiCorrected[iG]);

            if (options.umiDedup.yes.Directional)
                countCellGeneUMI[countCellGeneUMIindex[iCB + 1] +
                                 options.umiDedup.countInd.Directional] =
                    correctUmisDirectional(nU0, umiArray.data(),
                                           readInfo.size() > 0 &&
                                               options.umiDedup.typeMain ==
                                                   UmiDedupOptions::typeI::Directional,
                                           true, umiCorrected[iG], 0);

            if (options.umiDedup.yes.Directional_UMItools)
                countCellGeneUMI[countCellGeneUMIindex[iCB + 1] +
                                 options.umiDedup.countInd.Directional_UMItools] =
                    correctUmisDirectional(nU0, umiArray.data(),
                                           readInfo.size() > 0 &&
                                               options.umiDedup.typeMain ==
                                                   UmiDedupOptions::typeI::Directional_UMItools,
                                           true, umiCorrected[iG], -1);

            // this changes umiArray, so it should be last call
            if (options.umiDedup.yes.All)
                countCellGeneUMI[countCellGeneUMIindex[iCB + 1] + options.umiDedup.countInd.All] =
                    correctUmisGraph(nU0, umiArray.data(),
                                     readInfo.size() > 0 &&
                                         options.umiDedup.typeMain == UmiDedupOptions::typeI::All,
                                     true, umiCorrected[iG]);
        }; // if (nU0>0)

        { // check any count>0 and finalize record for this gene
            std::uint32_t totcount = 0;
            for (std::uint32_t ii = countCellGeneUMIindex[iCB + 1] + 1;
                 ii < countCellGeneUMIindex[iCB + 1] + countMatStride; ii++) {
                totcount += countCellGeneUMI[ii];
            };
            if (totcount > 0) { // at least one umiDedup type is non-0
                countCellGeneUMI[countCellGeneUMIindex[iCB + 1] + 0] = gID[iG];
                nGenePerCB[iCB]++;
                nUMIperCB[iCB] += countCellGeneUMI[countCellGeneUMIindex[iCB + 1] +
                                                   options.umiDedup.countInd.main];
                countCellGeneUMIindex[iCB + 1] =
                    countCellGeneUMIindex[iCB + 1] + countMatStride; // iCB+1 accumulates the index
            };
        };

        if (readInfo.size() > 0) { // record cb/umi for each read
            for (std::uint32_t iR = 0; iR < gReadS[iG + 1] - gReadS[iG];
                 iR += rguStride) { // cycle over reads
                UnsignedWord iread1 = rGU1[iR + rguR];
                readInfo[iread1].cb = indCB[iCB];
                std::uint32_t umi = rGU1[iR + rguU];

                if (umiCorrected[iG].count(umi) > 0)
                    umi = umiCorrected[iG][umi]; // correct UMI
                readInfo[iread1].umi = umi;
            };
        };
    };

    if (options.umiFiltering.MultiGeneUMI_CR) {

        std::vector<std::uint32_t> geneCounts(nGenes, 0);

        std::vector<std::unordered_set<UmiCode>> geneUmiHash;
        if (readInfo.size() > 0)
            geneUmiHash.resize(nGenes);

        for (const auto& iu : umiGeneMapCount) { // loop over UMIs for all genes

            std::uint32_t maxu = 0, maxg = -1;
            for (const auto& ig : iu.second) {
                if (ig.second > maxu) {
                    maxu = ig.second;
                    maxg = ig.first;
                } else if (ig.second == maxu) {
                    maxg = -1;
                };
            };

            if (maxg + 1 == 0)
                continue; // this umi is not counted for any gene, because two genes have the same
                          // read count for this UMI

            for (const auto& ig : umiGeneMapCount0[iu.first]) { // check that this umi/gene had also
                                                                // top count for uncorrected umis
                if (ig.second > umiGeneMapCount0[iu.first][maxg]) {
                    maxg = -1;
                    break;
                };
            };

            if (maxg + 1 != 0) { // this UMI is counted
                geneCounts[maxg]++;
                if (readInfo.size() > 0)
                    geneUmiHash[maxg].insert(iu.first);
            };
        };

        for (std::uint32_t ig = 0; ig < nGenes; ig++) {
            if (geneCounts[ig] == 0)
                continue; // no counts for this gene
            nGenePerCB[iCB]++;
            nUMIperCB[iCB] += geneCounts[ig];
            countCellGeneUMI[countCellGeneUMIindex[iCB + 1] + 0] = gID[ig];
            countCellGeneUMI[countCellGeneUMIindex[iCB + 1] + options.umiDedup.countInd.CR] =
                geneCounts[ig];
            countCellGeneUMIindex[iCB + 1] =
                countCellGeneUMIindex[iCB + 1] + countMatStride; // iCB+1 accumulates the index
        };

        if (readInfo.size() > 0) {                          // record cb/umi for each read
            for (std::uint32_t iG = 0; iG < nGenes; iG++) { // cycle over genes
                std::uint32_t* rGU1 = rGU + gReadS[iG];
                for (std::uint32_t iR = 0; iR < gReadS[iG + 1] - gReadS[iG];
                     iR += rguStride) { // cycle over reads
                    UnsignedWord iread1 = rGU1[iR + rguR];
                    readInfo[iread1].cb = indCB[iCB];
                    std::uint32_t umi = rGU1[iR + rguU];

                    if (umiCorrected[iG].count(umi) > 0)
                        umi = umiCorrected[iG][umi]; // correct UMI

                    // cout << iG << "-" << iR << " " <<flush ;
                    if (geneUmiHash[iG].count(umi) > 0) {
                        readInfo[iread1].umi = umi;
                    } else {
                        readInfo[iread1].umi = (UmiCode)-1;
                    };
                };
            };
        };
    };
    if (options.multiMap.yes.multi)
        countMatMult.i[iCB + 1] = countMatMult.i[iCB];

    if (nGenesMult > 0) { // process multigene reads

        if (readInfo.size() > 0) {
            for (std::uint32_t iR = gReadS[nGenes]; iR < gReadS[nGenes + nGenesMult];
                 iR += rguStride) { // cycle over multi-gene reads to record their CB and UMI, no
                                    // corrections
                UnsignedWord iread1 = rGU[iR + rguR];
                readInfo[iread1].cb = indCB[iCB];
                readInfo[iread1].umi = rGU[iR + rguU];
            };
        };

        std::vector<std::vector<std::uint32_t>> umiGenes;
        umiGenes.reserve(256);
        { // for each umi, count number of reads per gene.
          // Output umiGenes: only genes with nReads = nReads-for-this-UMI will be kept for this UMI
            std::uint32_t* rGUm = rGU + gReadS[nGenes];
            std::uint32_t nRm = (gReadS[nGenes + nGenesMult] - gReadS[nGenes]) / rguStride;

            // sort by UMI, then by read, then by gene
            qsort(rGUm, nRm, rguStride * sizeof(std::uint32_t),
                  funCompare_uint32_1_2_0); // there is no need to sort by read or gene actually

            std::unordered_map<std::uint32_t, std::uint32_t>
                geneReadCount; // number of reads per gene
            std::uint32_t nRumi = 0;
            bool skipUMI = false;
            UmiCode umiPrev = (UmiCode)-1;
            ReadCount readPrev = (ReadCount)-1;
            for (std::uint32_t iR = 0; iR < nRm * rguStride;
                 iR += rguStride) { // for each umi, find intersection of genes from each read
                UmiCode umi1 = rGUm[iR + 1];
                if (umi1 != umiPrev) { // starting new UMI
                    umiPrev = umi1;
                    if (umiGeneMapCount.count(umi1) > 0) {
                        skipUMI = true; // this UMI is skipped because it was among uniquely mapped
                    } else {
                        skipUMI = false; // new good umi
                        geneReadCount.clear();
                        nRumi = 0;
                        readPrev = (ReadCount)-1;
                    };
                };

                if (skipUMI)
                    continue; // this UMI is skipped because it was among uniquely mapped

                ReadCount read1 = rGUm[iR + 2];
                if (read1 != readPrev) {
                    ++nRumi;
                    readPrev = read1;
                };

                std::uint32_t g1 = rGUm[iR + 0] ^ geneMultMark; // XOR to unset the geneMultMark bit
                geneReadCount[g1]++;

                if (iR == nRm * rguStride - rguStride ||
                    umi1 != rGUm[iR + 1 + rguStride]) { // record this umi
                    std::uint32_t ng = 0;
                    for (const auto& gg : geneReadCount) {
                        if (gg.second == nRumi)
                            ++ng;
                    };
                    std::vector<std::uint32_t> vg;
                    vg.reserve(ng); // this and above is to construct vector of precise size, for
                                    // efficiency?
                    for (const auto& gg : geneReadCount) {
                        if (gg.second == nRumi)
                            vg.push_back(gg.first);
                    };
                    umiGenes.push_back(vg);
                };
            };
        };

        std::map<std::uint32_t, std::uint32_t> genesM; // genes to quantify

        { // collect all genes, replace geneID with index in umiGenes
            std::uint32_t ng = 0;
            for (auto& uu : umiGenes) {
                for (auto& gg : uu) {
                    if (genesM.count(gg) == 0) { // new gene
                        genesM[gg] = ng;
                        ++ng;
                    };
                    gg = genesM[gg];
                };
            };
        };

        std::vector<double> gEuniform(genesM.size(), 0);
        { // gEuniform=uniformly distribute multigene UMIs
            for (auto& ug : umiGenes) {
                for (auto& gg : ug) {
                    gEuniform[gg] += 1.0 / double(ug.size()); // 1/n_genes_umi
                };
            };
        };

        std::vector<std::vector<double>> gErescue(options.umiDedup.yes.N);
        if (options.multiMap.yes.Rescue) {

            for (std::uint32_t indDedup = 0; indDedup < options.umiDedup.yes.N; indDedup++) {
                std::vector<double> gEu(genesM.size(), 0);
                { // collect unique gene counts
                    for (std::uint32_t igm = countCellGeneUMIindex[iCB];
                         igm < countCellGeneUMIindex[iCB + 1]; igm += countMatStride) {
                        std::uint32_t g1 = countCellGeneUMI[igm];
                        if (genesM.count(g1) > 0)
                            gEu[genesM[g1]] = (double)countCellGeneUMI[igm + 1 + indDedup];
                    };
                };

                gErescue[indDedup].resize(genesM.size(), 0);
                { // gErescue=distribute UMI proportionally to gEuniform+gEu
                    for (auto& ug : umiGenes) {
                        double norm1 = 0.0;
                        for (auto& gg : ug)
                            norm1 += gEuniform[gg] + gEu[gg];

                        if (norm1 == 0.0)
                            continue; // this should not happen since gEuniform is non-zero for all
                                      // genes involved
                        norm1 = 1.0 / norm1;

                        for (auto& gg : ug) {
                            gErescue[indDedup][gg] += (gEuniform[gg] + gEu[gg]) * norm1;
                        };
                    };
                };
            };
        };

        std::vector<std::vector<double>> gEpropUnique(options.umiDedup.yes.N);
        if (options.multiMap.yes.PropUnique) {

            for (std::uint32_t indDedup = 0; indDedup < options.umiDedup.yes.N; indDedup++) {
                std::vector<double> gEu(genesM.size(), 0);
                { // collect unique gene counts
                    for (std::uint32_t igm = countCellGeneUMIindex[iCB];
                         igm < countCellGeneUMIindex[iCB + 1]; igm += countMatStride) {
                        std::uint32_t g1 = countCellGeneUMI[igm];
                        if (genesM.count(g1) > 0)
                            gEu[genesM[g1]] = (double)countCellGeneUMI[igm + 1 + indDedup];
                    };
                };

                gEpropUnique[indDedup].resize(genesM.size(), 0);
                { // gErescue=distribute UMI proportionally to gEuniform+gEu
                    for (auto& ug : umiGenes) {
                        double norm1 = 0.0;
                        for (auto& gg : ug)
                            norm1 += gEu[gg];

                        if (norm1 == 0.0) { // this UMI has no genes with unique mappers -
                                            // distribute it uniformly
                            for (auto& gg : ug)
                                gEpropUnique[indDedup][gg] += 1.0 / double(ug.size());
                        } else { // this UMI has genes with unique mappers - distribute it
                                 // proportionally to unique mappers
                            norm1 = 1.0 / norm1;
                            for (auto& gg : ug)
                                gEpropUnique[indDedup][gg] += gEu[gg] * norm1;
                        };
                    };
                };
            };
        };

        std::vector<std::vector<double>> gEem(options.umiDedup.yes.N);
        if (options.multiMap.yes.EM) {

            for (std::uint32_t indDedup = 0; indDedup < options.umiDedup.yes.N; indDedup++) {
                std::vector<double> gEu(genesM.size(), 0);
                { // collect unique gene counts
                    for (std::uint32_t igm = countCellGeneUMIindex[iCB];
                         igm < countCellGeneUMIindex[iCB + 1]; igm += countMatStride) {
                        std::uint32_t g1 = countCellGeneUMI[igm];
                        if (genesM.count(g1) > 0)
                            gEu[genesM[g1]] = (double)countCellGeneUMI[igm + 1 + indDedup];
                    };
                };

                { // gEem = EM

                    std::vector<double> gEM1 = gEuniform;
                    for (std::uint32_t ii = 0; ii < gEM1.size(); ii++)
                        gEM1[ii] += gEu[ii]; // start with sum of unique and uniform

                    std::vector<double> gEM2(genesM.size(), 0);

                    auto* gEM1p = &gEM1;
                    auto* gEM2p = &gEM2;

                    double maxAbsChange = 1;
                    std::uint32_t iterI = 0;
                    while (true) {
                        ++iterI;

                        auto& gEMold = *gEM1p; // for convenience - to use instead of pointer
                        auto& gEMnew = *gEM2p;

                        std::copy(gEu.begin(), gEu.end(),
                                  gEMnew.begin()); // gEMnew is initialized with unique counts

                        for (auto& gg : gEMold) { // zero-out very small counts
                            if (gg < 0.01)        // hardcoded
                                gg = 0;
                        };

                        for (auto& ug : umiGenes) {
                            double norm1 = 0.0;
                            for (auto& gg : ug) // cycle over genes for this umi
                                norm1 += gEMold[gg];

                            norm1 = 1.0 / norm1;

                            for (auto& gg : ug) {
                                gEMnew[gg] += gEMold[gg] * norm1;
                            };
                        };

                        maxAbsChange = 0.0;
                        for (std::uint32_t ii = 0; ii < gEMnew.size(); ii++) {
                            double change1 = abs(gEMnew[ii] - gEMold[ii]);
                            if (change1 > maxAbsChange)
                                maxAbsChange = change1;
                        };

                        if (maxAbsChange < 0.01 || iterI > 100) { // hardcoded
                            gEem[indDedup] = gEMnew;
                            break;
                        };

                        std::swap(gEM1p, gEM2p); // swap new and old for the next interation of EM
                    };

                    for (std::uint32_t ii = 0; ii < gEM1.size(); ii++)
                        gEem[indDedup][ii] -= gEu[ii]; // gEem contains only multimapper counts
                };
            };
        };

        { // write to countMatMult
            if (countMatMult.m.size() <
                countMatMult.i[iCB + 1] + genesM.size() * countMatMult.s * options.umiDedup.yes.N +
                    100) //+100 just in case
                countMatMult.m.resize((countMatMult.i[iCB + 1] +
                                       genesM.size() * countMatMult.s * options.umiDedup.yes.N +
                                       100) *
                                      2);

            for (const auto& gm : genesM) {
                countMatMult.m[countMatMult.i[iCB + 1] + 0] = gm.first;

                for (std::uint32_t indDedup = 0; indDedup < options.umiDedup.yes.N; indDedup++) {
                    std::uint32_t ind1 = countMatMult.i[iCB + 1] + indDedup;

                    if (options.multiMap.yes.Uniform)
                        countMatMult.m[ind1 + options.multiMap.countInd.Uniform] =
                            gEuniform[gm.second];

                    if (options.multiMap.yes.Rescue)
                        countMatMult.m[ind1 + options.multiMap.countInd.Rescue] =
                            gErescue[indDedup][gm.second];

                    if (options.multiMap.yes.PropUnique)
                        countMatMult.m[ind1 + options.multiMap.countInd.PropUnique] =
                            gEpropUnique[indDedup][gm.second];

                    if (options.multiMap.yes.EM)
                        countMatMult.m[ind1 + options.multiMap.countInd.EM] =
                            gEem[indDedup][gm.second];

                    countMatMult.i[iCB + 1] += countMatMult.s;
                };
            };
        };
    };
};
inline int funCompareSolo1(const void* a, const void* b) {
    std::uint32_t* va = (std::uint32_t*)a;
    std::uint32_t* vb = (std::uint32_t*)b;

    if (va[1] > vb[1]) {
        return 1;
    } else if (va[1] < vb[1]) {
        return -1;
    } else if (va[0] > vb[0]) {
        return 1;
    } else if (va[0] < vb[0]) {
        return -1;
    } else {
        return 0;
    };
};

inline int funCompare_uint32_1_2_0(const void* a, const void* b) {
    std::uint32_t* va = (std::uint32_t*)a;
    std::uint32_t* vb = (std::uint32_t*)b;

    if (va[1] > vb[1]) {
        return 1;
    } else if (va[1] < vb[1]) {
        return -1;
    } else if (va[2] > vb[2]) {
        return 1;
    } else if (va[2] < vb[2]) {
        return -1;
    } else if (va[0] > vb[0]) {
        return 1;
    } else if (va[0] < vb[0]) {
        return -1;
    } else {
        return 0;
    };
};
std::uint32_t FeatureCounter::correctUmisCellRanger(const std::uint32_t nU0, UmiCode* umiArr,
                                                    const bool readInfoRec, const bool nUMIyes,
                                                    std::unordered_map<UmiCode, UmiCode>& umiCorr) {
    qsort(umiArr, nU0, umiArrayStride * sizeof(std::uint32_t), funCompareSolo1);

    for (UnsignedWord iu = 0; iu < nU0 * umiArrayStride; iu += umiArrayStride) {

        umiArr[iu + 2] = umiArr[iu + 0]; // stores corrected UMI for 1MM_CR and 1MM_Directional
        for (UnsignedWord iuu = (nU0 - 1) * umiArrayStride; iuu > iu; iuu -= umiArrayStride) {

            std::uint32_t uuXor = umiArr[iu + 0] ^ umiArr[iuu + 0];

            if ((uuXor >> (__builtin_ctz(uuXor) / 2) * 2) <= 3) { // 1MM
                umiArr[iu + 2] = umiArr[iuu + 0];                 // replace iu with iuu
                break;
            };
        };
    };

    if (readInfoRec) { // record corrections
        for (UnsignedWord iu = 0; iu < nU0 * umiArrayStride; iu += umiArrayStride) {
            if (umiArr[iu + 0] != umiArr[iu + 2])
                umiCorr[umiArr[iu + 0]] = umiArr[iu + 2];
        };
    };

    if (!nUMIyes) {
        return 0;
    } else {
        std::unordered_set<UmiCode> umiC;
        for (UnsignedWord iu = 0; iu < nU0 * umiArrayStride; iu += umiArrayStride) {
            umiC.insert(umiArr[iu + 2]);
        };
        return umiC.size();
    };
};
std::uint32_t FeatureCounter::correctUmisDirectional(const std::uint32_t nU0, UmiCode* umiArr,
                                                     const bool readInfoRec, const bool nUMIyes,
                                                     std::unordered_map<UmiCode, UmiCode>& umiCorr,
                                                     const std::int32_t dirCountAdd) {
    qsort(umiArr, nU0, umiArrayStride * sizeof(std::uint32_t),
          funCompareNumbersReverseShift<std::uint32_t, 1>); // TODO no need to sort by sequence
                                                            // here, only by count.

    for (UnsignedWord iu = 0; iu < nU0 * umiArrayStride; iu += umiArrayStride)
        umiArr[iu + 2] =
            umiArr[iu +
                   0]; // initialized - it will store corrected UMI for 1MM_CR and 1MM_Directional

    std::uint32_t nU1 = nU0;
    for (UnsignedWord iu = umiArrayStride; iu < nU0 * umiArrayStride; iu += umiArrayStride) {

        for (UnsignedWord iuu = 0; iuu < iu; iuu += umiArrayStride) {

            std::uint32_t uuXor = umiArr[iu + 0] ^ umiArr[iuu + 0];

            if ((uuXor >> (__builtin_ctz(uuXor) / 2) * 2) <= 3 &&
                umiArr[iuu + 1] >=
                    (2 * umiArr[iu + 1] + dirCountAdd)) { // 1MM && directional condition
                umiArr[iu + 2] = umiArr[iuu + 2];         // replace iuu with iu-corrected
                nU1--;
                break;
            };
        };
    };

    if (readInfoRec) { // record corrections
        for (UnsignedWord iu = 0; iu < nU0 * umiArrayStride; iu += umiArrayStride) {
            if (umiArr[iu + 0] != umiArr[iu + 2])
                umiCorr[umiArr[iu + 0]] = umiArr[iu + 2];
        };
    };

    if (!nUMIyes) {
        return 0;
    } else {
        std::unordered_set<UmiCode> umiC;
        for (UnsignedWord iu = 0; iu < nU0 * umiArrayStride; iu += umiArrayStride) {
            umiC.insert(umiArr[iu + 2]);
        };
        if (umiC.size() != nU1)
            std::cout << nU1 << " " << umiC.size() << std::endl;
        return umiC.size();
    };
};

} // namespace nustar::star2711b::solo::detail
