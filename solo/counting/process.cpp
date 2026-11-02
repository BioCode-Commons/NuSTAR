#include "../barcode/sequence.hpp"
#include "feature_counter.hpp"
#include "../support/error.hpp"
#include "../support/memory.hpp"
#include "../support/streams.hpp"
#include "../support/time.hpp"

namespace nustar::star2711b::solo::detail {

void FeatureCounter::processRecords() {
    if (options.type == 0)
        return;

    time_t rawTime;
    time(&rawTime);
    context_.inOut->logMain << timeMonthDayTime(rawTime) << " ... Starting Solo post-map for "
                            << FeatureTypes::Names[featureType] << std::endl;

    outputPrefix = context_.outFileNamePrefix + options.outFileNames[0];
    outputPrefix += FeatureTypes::Names[featureType] + '/';
    outputPrefixFiltered = outputPrefix + "filtered/";

    if (mkdir(outputPrefix.c_str(), context_.runDirPerm) != 0 &&
        errno != EEXIST) { // create directory
        std::ostringstream errOut;
        errOut
            << "EXITING because of fatal OUTPUT FILE error: could not create Solo output directory"
            << outputPrefix << "\n";
        errOut << "SOLUTION: check the path and permisssions";
        exitWithError(errOut.str(), std::cerr, context_.inOut->logMain, EXIT_CODE_PARAMETER,
                      context_);
    };

    // prepare for feature-specific counting:
    if (featureType == FeatureTypes::SJ && context_.sjAll[0].empty()) {
        std::ifstream& sjIn = ifstrOpen(context_.outFileTmp + "SJ.start_gap.tsv", sourceLocation(),
                                        "SOLUTION: re-run STAR", context_);
        context_.sjAll[0].reserve(10000000);
        context_.sjAll[1].reserve(10000000);
        UnsignedWord start1, gap1;
        while (sjIn >> start1 >> gap1) {
            context_.sjAll[0].emplace_back(start1);
            context_.sjAll[1].emplace_back(gap1);
        };
        sjIn.close();
        context_.inOut->logMain << "Read splice junctions for Solo SJ feature: "
                                << context_.sjAll[0].size() << std::endl;
    };

    FeatureCounter::sumThreads();

    // call counting method
    if (featureType == FeatureTypes::Velocyto) {
        countVelocyto();
    } else if (featureType == FeatureTypes::Transcript3p) {
        quantTranscript();
        return;
    } else { // all others, standard processing
        if (options.type == options.SoloTypes::SmartSeq) {
            countSmartSeq();
        } else {
            countCBgeneUMI();
        };
    };

    // no need to include multi-gene in yesWLmatch:
    // for multi-gene options, it's already done in _collapseUMI_ function
    //    readFeatSum->stats.V[readFeatSum->stats.yesWLmatch] +=
    //    readFeatSum->stats.V[readFeatSum->stats.MultiFeature];

    // output
    std::ofstream* statsStream =
        &ofstrOpen(outputPrefix + "Features.stats", sourceLocation(), context_);
    readFeatSum->statsOut(*statsStream);
    statsStream->close();

    time(&rawTime);
    context_.inOut->logMain << timeMonthDayTime(rawTime) << " ... Solo: writing raw matrix"
                            << std::endl;

    // output nU per gene per CB
    outputResults(false, outputPrefix + "/raw/"); // unfiltered

    time(&rawTime);
    context_.inOut->logMain << timeMonthDayTime(rawTime) << " ... Solo: cell filtering"
                            << std::endl;
    cellFiltering();

    // summary stats output
    statsOutput();
    clearLarge();

    context_.inOut->logMain << "RAM after completing solo:\n" << linuxProcMemory() << std::flush;
};

} // namespace nustar::star2711b::solo::detail
