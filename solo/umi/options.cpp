#include "../configuration/options.hpp"
#include "../barcode/sequence.hpp"
#include "../configuration/context.hpp"
#include "../support/algorithms.hpp"
#include "../support/error.hpp"
#include "../support/streams.hpp"
#include <stdlib.h>

namespace nustar::star2711b::solo::detail {
void Options::umiSwapHalves(std::uint32_t& umi) {
    std::uint32_t high = umi >> (umiL);
    umi &= umiMaskLow; // remove high
    umi <<= (umiL);    // move low to high
    umi |= high;       // add high
};

void UmiDedupOptions::initialize(Options* pS) {
    yes.N = 0;
    countInd.I.fill((uint32_t)-1); // marks types not used
    yes.B.fill(false);

    for (uint32_t iin = 0; iin < typesIn.size(); iin++) {
        uint32_t itype;
        for (itype = 0; itype < tN; itype++) {
            if (typesIn[iin] == typeNames[itype])
                break; // found match
        };

        if (itype == tN) { // no match
            std::string tall;
            for (auto& t : typeNames)
                tall += " " + t; // concatenate allowed values

            exitWithError(
                "EXITING because of fatal PARAMETERS error: unrecognzied option --soloUMIdedup = " +
                    typesIn[iin] + '\n' + "SOLUTION: use allowed values: " + tall + '\n',
                std::cerr, pS->context->inOut->logMain, EXIT_CODE_PARAMETER, *pS->context);
        };

        types.push_back(itype);
        yes.B[itype] = true;
        yes.N++;
        countInd.I[itype] = iin + 1; // for each type, which column itype's recorded in

        if (pS->type == pS->SoloTypes::SmartSeq && (yes.All || yes.Directional || yes.CR))
            exitWithError(
                "EXITING because of fatal PARAMETERS error: --soloUMIdedup = " + typesIn[iin] +
                    " is not allowed for --soloType SmartSeq\n" +
                    "SOLUTION: use allowed options: Exact and/or NoDedup\n",
                std::cerr, pS->context->inOut->logMain, EXIT_CODE_PARAMETER, *pS->context);
    };

    // hard-coded for now
    typeMain = types[0]; // main is the 0th entry in typesIn
    countInd.main = 1;   // hard-coded - 1 column is always main
};

} // namespace nustar::star2711b::solo::detail
