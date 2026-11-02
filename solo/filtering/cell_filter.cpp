#include "../counting/feature_counter.hpp"
#include "../support/algorithms.hpp"
#include <math.h>

namespace nustar::star2711b::solo::detail {

void FeatureCounter::cellFiltering() {

    if (options.cellFilter.type[0] == "None" || nCB < 1) { // no filtering, or no cells to filter
        return;
    };

    switch (featureType) {
    default:
        return; // filtering not done for other features

    case FeatureTypes::Velocyto: {
        filteredCells.reset(nCB);

        // Velocyto: use filter from Gene
        FeatureCounter& soFeGe = *features_all_[options.featureInd[FeatureTypes::Gene]];
        for (std::uint32_t ic = 0; ic < soFeGe.nCB; ic++) {
            if (soFeGe.filteredCells.filtVecBool[ic] &&
                indCBwl[soFeGe.indCB[ic]] != (std::uint32_t)-1) {
                // this cell passde Gene filtering, and is detected in Velocyto
                filteredCells.filtVecBool[indCBwl[soFeGe.indCB[ic]]] = true;
            };
        };

        nUMIperCBsorted = nUMIperCB;
        std::sort(
            nUMIperCBsorted.begin(), nUMIperCBsorted.end(),
            [](const uint32_t& u1, const uint32_t& u2) { return u1 > u2; }); // sort descending
        break;
    };
    case FeatureTypes::Gene:
    case FeatureTypes::GeneFull:
    case FeatureTypes::GeneFull_Ex50pAS:
    case FeatureTypes::GeneFull_ExonOverIntron:
    case -1: // undefined: matrix loaded from file
    {
        // do the filtering

        // simple filtering first
        nUMIperCBsorted = nUMIperCB;
        std::sort(
            nUMIperCBsorted.begin(), nUMIperCBsorted.end(),
            [](const uint32_t& u1, const uint32_t& u2) { return u1 > u2; }); // sort descending

        std::uint32_t nUMImax = 0, nUMImin = 0;
        if (options.cellFilter.type[0] == "TopCells") {
            nUMImin = nUMIperCBsorted[std::min(nCB - 1, options.cellFilter.topCells)];
        } else { // other filtering types require simple filtering first
            // find robust max
            std::uint32_t maxind = int(
                round(options.cellFilter.knee.nExpectedCells *
                      (1.0 - options.cellFilter.knee.maxPercentile))); // cell number for robust max
            nUMImax = nUMIperCBsorted[std::min(nCB - 1, maxind)]; // robust estimate of the max UMI
            nUMImin = int(round(nUMImax / options.cellFilter.knee.maxMinRatio));
        };
        nUMImin = std::max(nUMImin, (std::uint32_t)1); // cannot be zero

        filteredCells.reset(nCB); // all stats to 0

        for (std::uint32_t icb = 0; icb < nCB; icb++) {
            if (nUMIperCB[icb] >= nUMImin) {
                filteredCells.filtVecBool[icb] = true;
                filteredCells.nCellsSimple++;
            };
        };

        context_.inOut->logMain << "cellFiltering: simple: nUMImax=" << nUMImax
                                << "; nUMImin=" << nUMImin
                                << "; nCellsSimple=" << filteredCells.nCellsSimple << std::endl;

        if (options.cellFilter.type[0] == "EmptyDrops_CR") {
            emptyDrops_CR();
        };
    };
    };
    // filtering is done: filtVecBool=true for kept cells

    // calculate filtered statistics

    std::vector<uint32_t> geneDetected(featuresNumber,
                                       0); //=1 if a gene was detected in at least one cell

    for (std::uint32_t icb = 0; icb < nCB; icb++) {
        if (filteredCells.filtVecBool[icb]) {

            filteredCells.nCells++;

            filteredCells.nUMIinCells +=
                nUMIperCB[icb]; // nUMIperCB was calculated for umiDedup-main

            if (nReadPerCBunique.size() >
                0) { // for CellFiltering only, read information is not available
                filteredCells.nReadInCells += nReadPerCBtotal[icb];
                filteredCells.nReadInCellsUnique += nReadPerCBunique[icb];
                filteredCells.nReadPerCellUnique.push_back(nReadPerCBunique[icb]);
            };

            std::uint32_t ng1 = 0; // number of genes detected in this cell
            for (std::uint32_t ig = 0; ig < nGenePerCB[icb]; ig++) {
                std::uint32_t indG1 = countCellGeneUMIindex[icb] + ig * countMatStride;
                if (countCellGeneUMI[indG1 + options.umiDedup.countInd.main] > 0) {
                    geneDetected[countCellGeneUMI[indG1]] =
                        1; // gene is present if it's count > 0 for
                    ng1++;
                };
            };

            filteredCells.nGeneInCells +=
                ng1; // had to recalculate this since some gene counts could be 0
            filteredCells.nGenePerCell.push_back(ng1);
        };
    };

    if (filteredCells.nCells == 0) { // all stats were already set to 0
        return;
    };

    filteredCells.nGeneDetected = 0;
    for (std::uint32_t ii = 0; ii < featuresNumber; ii++) {
        if (geneDetected[ii] > 0)
            filteredCells.nGeneDetected++;
    };

    filteredCells.meanUMIperCell = filteredCells.nUMIinCells / filteredCells.nCells;
    filteredCells.meanReadPerCellUnique = filteredCells.nReadInCellsUnique / filteredCells.nCells;
    filteredCells.meanGenePerCell = filteredCells.nGeneInCells / filteredCells.nCells;

    std::sort(filteredCells.nReadPerCellUnique.begin(), filteredCells.nReadPerCellUnique.end());
    std::sort(filteredCells.nGenePerCell.begin(), filteredCells.nGenePerCell.end());

    filteredCells.medianUMIperCell = nUMIperCBsorted[filteredCells.nCells / 2];
    filteredCells.medianGenePerCell = filteredCells.nGenePerCell[filteredCells.nCells / 2];
    // Matrix-only filtering has no per-cell read counts. Upstream indexes an
    // empty vector here; this statistic is not written in matrix-only mode.
    if (context_.runMode != "soloCellFiltering")
        filteredCells.medianReadPerCellUnique =
            filteredCells.nReadPerCellUnique[filteredCells.nCells / 2];
    outputResults(true, outputPrefixFiltered);

    return;
};

} // namespace nustar::star2711b::solo::detail
