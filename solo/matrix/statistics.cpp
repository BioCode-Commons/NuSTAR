#include "../counting/feature_counter.hpp"
#include "../support/streams.hpp"

namespace nustar::star2711b::solo::detail {

void FeatureCounter::statsOutput() {
    std::ofstream& strOut = ofstrOpen(outputPrefix + "Summary.csv", sourceLocation(), context_);
    // Sequencing
    strOut << "Number of Reads," << context_.countingStats.readN << '\n';
    strOut << "Reads With Valid Barcodes,"
           << 1.0 - double(barcode_sum_->stats.numInvalidBarcodes() +
                           readFeatSum->stats.numInvalidBarcodes()) /
                        context_.countingStats.readN
           << '\n';
    strOut << "Sequencing Saturation," << readFeatSum->stats.numSequencingSaturation() << '\n';

    if (options.type != options.SoloTypes::SmartSeq) { // qualHist for CB+UMI
        UnsignedWord q30 = 0, ntot = 0;
        for (std::uint32_t ix = 0; ix < 256; ix++) {
            ntot += barcode_sum_->qualHist[ix];
            if (ix >= (context_.readQualityScoreBase + 30))
                q30 += barcode_sum_->qualHist[ix];
        };

        strOut << "Q30 Bases in CB+UMI," << double(q30) / ntot << '\n';
    };

    { // qualHist for RNA
        UnsignedWord q30 = 0, ntot = 0;
        for (int ichunk = 0; ichunk < context_.runThreadN; ichunk++) {
            for (std::uint32_t imate = 0; imate < context_.readNmates; imate++) {
                for (std::uint32_t ix = 0; ix < 256; ix++) {
                    ntot += workers[ichunk]->qualHist[imate][ix];
                    if (ix >= (context_.readQualityScoreBase + 30))
                        q30 += workers[ichunk]->qualHist[imate][ix];
                };
            };
        };
        strOut << "Q30 Bases in RNA read," << double(q30) / ntot << '\n';
    };

    strOut << "Reads Mapped to Genome: Unique+Multiple,"
           << double(context_.countingStats.mappedReadsU + context_.countingStats.mappedReadsM) /
                  context_.countingStats.readN
           << '\n';
    strOut << "Reads Mapped to Genome: Unique,"
           << double(context_.countingStats.mappedReadsU) / context_.countingStats.readN << '\n';

    std::string mapfeat = FeatureTypes::Names[featureType];

    strOut << "Reads Mapped to " << mapfeat << ": Unique+Multiple "
           << FeatureTypes::Names[featureType] << ",";
    if (options.multiMap.yes.multi) {
        strOut << double(readFeatSum->stats.numMappedToTranscriptome()) /
                      context_.countingStats.readN
               << '\n';
    } else {
        strOut << "NoMulti\n";
    };

    strOut << "Reads Mapped to " << mapfeat << ": Unique " << FeatureTypes::Names[featureType]
           << ","
           << double(readFeatSum->stats.numMappedToTranscriptomeUnique()) /
                  context_.countingStats.readN
           << '\n';

    if (options.cellFilter.type[0] != "None" &&
        (featureType == FeatureTypes::Gene || featureType == FeatureTypes::GeneFull ||
         featureType == FeatureTypes::GeneFull_ExonOverIntron ||
         featureType == FeatureTypes::GeneFull_Ex50pAS)) {
        {
            strOut << "Estimated Number of Cells," << filteredCells.nCells << '\n';

            strOut << "Unique Reads in Cells Mapped to " << FeatureTypes::Names[featureType] << ","
                   << filteredCells.nReadInCellsUnique << '\n';
            strOut << "Fraction of Unique Reads in Cells,"
                   << double(filteredCells.nReadInCellsUnique) /
                          readFeatSum->stats.numMappedToTranscriptomeUnique()
                   << '\n';
            strOut << "Mean Reads per Cell," << filteredCells.meanReadPerCellUnique << '\n';
            strOut << "Median Reads per Cell," << filteredCells.medianReadPerCellUnique << '\n';

            strOut << "UMIs in Cells," << filteredCells.nUMIinCells << '\n';
            strOut << "Mean UMI per Cell," << filteredCells.meanUMIperCell << '\n';
            strOut << "Median UMI per Cell," << filteredCells.medianUMIperCell << '\n';

            strOut << "Mean " << FeatureTypes::Names[featureType] << " per Cell,"
                   << filteredCells.meanGenePerCell << '\n';
            strOut << "Median " << FeatureTypes::Names[featureType] << " per Cell,"
                   << filteredCells.medianGenePerCell << '\n';
            strOut << "Total " << FeatureTypes::Names[featureType] << " Detected,"
                   << filteredCells.nGeneDetected << '\n';
        };

        // output UMI per cell, sorted
        std::ofstream& strOutUMIperCell =
            ofstrOpen(outputPrefix + "UMIperCellSorted.txt", sourceLocation(), context_);

        for (auto& n : nUMIperCBsorted) {
            if (n == 0)
                break;
            strOutUMIperCell << n << '\n';
        };
        strOutUMIperCell.close();
    };

    strOut.close();
    if (options.readStatsYes[featureType]) {
        std::ofstream& strOut =
            ofstrOpen(outputPrefix + "CellReads.stats", sourceLocation(), context_);

        strOut << "CB";
        for (auto& sn : readFlagCounts.statNames)
            strOut << "\t" << sn;
        strOut << "\t" << "nUMIunique" << "\t" << "nGenesUnique" << "\t" << "nUMImulti" << "\t"
               << "nGenesMulti";
        strOut << '\n';

        strOut << "CBnotInPasslist";
        for (auto& cc : readFlagCounts.flagCountsNoCB)
            strOut << "\t" << cc;
        strOut << "\t0\t0\t0\t0\n";

        for (auto& cbc : readFlagCounts.flagCounts) {
            strOut << options.cbWLstr[cbc.first];

            for (std::uint32_t ib = 0; ib < readFlagCounts.nBits; ib++)
                strOut << '\t' << cbc.second[ib];

            if (indCBwl[cbc.first] ==
                (std::uint32_t)-1) { // this CB did not have any reads mapped to features
                strOut << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0;
            } else {
                strOut << '\t' << nUMIperCB[indCBwl[cbc.first]] << '\t'
                       << nGenePerCB[indCBwl[cbc.first]];
                if (nUMIperCBmulti.size() == 0) { // no multi
                    strOut << '\t' << 0 << '\t' << 0;
                } else {
                    strOut << '\t' << nUMIperCBmulti[indCBwl[cbc.first]] << '\t'
                           << nGenePerCBmulti[indCBwl[cbc.first]];
                };
            };
            strOut << '\n';
        };
        strOut.close();
    };
};

} // namespace nustar::star2711b::solo::detail
