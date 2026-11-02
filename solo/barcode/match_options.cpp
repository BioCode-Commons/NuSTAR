#include "sequence.hpp"
#include "../configuration/context.hpp"
#include "../configuration/options.hpp"
#include "../support/algorithms.hpp"
#include "../support/error.hpp"
#include "../support/streams.hpp"
#include <stdlib.h>

namespace nustar::star2711b::solo::detail {
void Options::init_CBmatchWL() { // CBmatchWL
    bool incomp1 =
        typeStr == "CB_UMI_Complex" &&
        (CBmatchWL.type != "Exact" && CBmatchWL.type != "1MM" && CBmatchWL.type != "EditDist_2");
    incomp1 = incomp1 ||
              (typeStr == "CB_samTagOut" && (CBmatchWL.type != "Exact" && CBmatchWL.type != "1MM"));
    incomp1 = incomp1 || (typeStr != "CB_UMI_Complex" && CBmatchWL.type == "EditDist_2");

    if (incomp1) {
        std::ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETERS error: --soloCBmatchWLtype "
               << CBmatchWL.type << " does not work with --soloType " << typeStr << "\n";
        errOut << "SOLUTION: use allowed option: use --soloCBmatchWLtype Exact (exact matches "
                  "only) OR 1MM (one match with 1 mismatched base)\n";
        exitWithError(errOut.str(), std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER,
                      *context);
    };

    CBmatchWL.mm1 = false;
    CBmatchWL.mm1_multi = false;
    CBmatchWL.mm1_multi_pc = false;
    CBmatchWL.mm1_multi_Nbase = false;
    CBmatchWL.oneExact =
        false; // if true, for a CB matching with 1 mismatch to a WL-CB, requires at least one other
               // read to match this WL-CB. this is true for all options except pseudocount
    CBmatchWL.EditDist_2 = false;

    if (CBmatchWL.type == "Exact") {
        CBmatchWL.oneExact = true;
    } else if (CBmatchWL.type == "1MM") {
        CBmatchWL.mm1 = true;
        CBmatchWL.oneExact = true;
    } else if (CBmatchWL.type == "1MM_multi") {
        CBmatchWL.mm1 = true;
        CBmatchWL.mm1_multi = true;
        CBmatchWL.oneExact = true;
    } else if (CBmatchWL.type == "1MM_multi_pseudocounts") {
        CBmatchWL.mm1 = true;
        CBmatchWL.mm1_multi = true;
        CBmatchWL.mm1_multi_pc = true;
    } else if (CBmatchWL.type == "1MM_multi_Nbase_pseudocounts") {
        CBmatchWL.mm1 = true;
        CBmatchWL.mm1_multi = true;
        CBmatchWL.mm1_multi_pc = true;
        CBmatchWL.mm1_multi_Nbase = true;
    } else if (CBmatchWL.type == "EditDist_2") {
        CBmatchWL.EditDist_2 = true;
    } else {
        exitWithError("EXITING because of fatal PARAMETERS error: unrecognized option in "
                      "--soloCBmatchWLtype " +
                          CBmatchWL.type +
                          "\nSOLUTION: use allowed options: Exact or 1MM or 1MM_multi or "
                          "1MM_multi_pseudocounts 1MM_multi_Nbase_pseudocounts\n",
                      std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER, *context);
    };
};

} // namespace nustar::star2711b::solo::detail
