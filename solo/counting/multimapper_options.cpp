#include "../barcode/sequence.hpp"
#include "../configuration/context.hpp"
#include "../configuration/options.hpp"
#include "../support/algorithms.hpp"
#include "../support/error.hpp"
#include "../support/streams.hpp"
#include <stdlib.h>

namespace nustar::star2711b::solo::detail {
void MultimapperOptions::initialize(Options* pS) {
    yes.N = 0;
    countInd.I.fill((std::uint32_t)-1); // marks types not used
    yes.B.fill(false);

    for (std::uint32_t iin = 0; iin < typesIn.size(); iin++) {
        std::uint32_t itype;
        for (itype = 0; itype < tN; itype++) {
            if (typesIn[iin] == typeNames[itype])
                break; // found match
        };

        if (itype == tN) { // no match
            std::string tall;
            for (auto& t : typeNames)
                tall += " " + t; // concatenate allowed values

            exitWithError("EXITING because of fatal PARAMETERS error: unrecognzied option "
                          "--soloMultiMappers = " +
                              typesIn[iin] + '\n' + "SOLUTION: use allowed values: " + tall + '\n',
                          std::cerr, pS->context->inOut->logMain, EXIT_CODE_PARAMETER,
                          *pS->context);
        };

        if (itype == typeI::Unique)
            continue; // Unique type does not have to be recorded TODO use it to perform filtering
                      // and stats

        types.push_back(itype);
        yes.B[itype] = true;
        yes.N++;
    };

    if (yes.N == 0) { // only Unique, no multimappers
        yes.multi = false;
        return;
    };

    std::uint32_t ind1 = 1; // start
    for (const auto& itype : types) {
        countInd.I[itype] = ind1;
        ind1 += pS->umiDedup.yes.N;
    };

    // hard-coded for now
    typeMain = types[0]; // main is the 0th entry in typesIn
    countInd.main = 1;   // hard-coded - 1 column is always main

    yes.multi = yes.Uniform | yes.Rescue | yes.PropUnique | yes.EM;
};

} // namespace nustar::star2711b::solo::detail
