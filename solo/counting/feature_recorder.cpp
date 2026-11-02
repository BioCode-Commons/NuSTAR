#include "feature_recorder.hpp"
#include "../annotation/feature_types.hpp"
#include "../support/streams.hpp"

namespace nustar::star2711b::solo::detail {

FeatureRecorder::FeatureRecorder(std::int32_t feTy, Context& context_in, int worker_id)
    : featureType(feTy), context_(context_in), options(context_.options) {
    if (options.type == 0)
        return;
    //         return;

    readInfoYes = options.readInfoYes[featureType];
    readIndexYes = options.readIndexYes[featureType];

    if (options.cbWLyes) {
        cbReadCount.resize(options.cbWLsize, 0);
    };

    if (worker_id >= 0) {
        // open with flagDelete=false, i.e. try to keep file if it exists
        streamReads = &fstrOpen(context_.outFileTmp + "/solo" + FeatureTypes::Names[featureType] +
                                    '_' + std::to_string(worker_id),
                                sourceLocation(), context_, false);
    };

    if (featureType == FeatureTypes::Transcript3p)
        transcriptDistCount.resize(10000, 0);
};

void FeatureRecorder::addCounts(const FeatureRecorder& rfIn) {
    if (options.cbWLyes) { // WL
        for (std::uint32_t ii = 0; ii < options.cbWLsize; ii++) {
            cbReadCount[ii] += rfIn.cbReadCount[ii];
        };
    } else {
        for (auto ii = rfIn.cbReadCountMap.cbegin(); ii != rfIn.cbReadCountMap.cend(); ++ii) {
            cbReadCountMap[ii->first] += ii->second;
        };
    };

    if (transcriptDistCount.size() > 0) {
        for (std::uint32_t ii = 0; ii < transcriptDistCount.size(); ii++)
            transcriptDistCount[ii] += rfIn.transcriptDistCount[ii];
    };
};

void FeatureRecorder::addStats(const FeatureRecorder& rfIn) {
    for (std::uint32_t ii = 0; ii < stats.nStats; ii++)
        stats.V[ii] += rfIn.stats.V[ii];

    for (std::uint32_t ii = 0; ii < readFlag.nBits; ii++)
        readFlag.flagCountsNoCB[ii] += rfIn.readFlag.flagCountsNoCB[ii];
};

void FeatureRecorder::statsOut(std::ofstream& streamOut) {
    // streamOut << setw(50) << "CELL BARCODES IN READS:\n"
    for (std::uint32_t ii = 0; ii < stats.nStats; ii++) {
        streamOut << std::setw(50) << stats.names[ii] << std::setw(15) << stats.V[ii] << '\n';
    };
    streamOut.flush();
};

} // namespace nustar::star2711b::solo::detail
