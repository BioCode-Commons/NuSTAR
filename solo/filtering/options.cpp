#include "../configuration/options.hpp"
#include "../barcode/sequence.hpp"
#include "../configuration/context.hpp"
#include "../support/algorithms.hpp"
#include "../support/error.hpp"
#include "../support/streams.hpp"
#include <stdlib.h>

namespace nustar::star2711b::solo::detail {
void Options::cellFiltering() { // cell filtering
    std::string pars1;
    for (auto s = cellFilter.type.begin() + 1; s != cellFilter.type.end(); s++)
        pars1 += ' ' + *s; // concatenate parameters into one string - easier to process that way

    if (cellFilter.type[0] == "CellRanger2.2") {
        if (cellFilter.type.size() == 1) {
            context->inOut->logMain << "ParametersSolo: using hardcoded filtering parameters for "
                                       "--soloCellFilter CellRanger2.2"
                                    << std::endl;
            pars1 = "3000 0.99 10";
        } else if (cellFilter.type.size() < 4) {
            std::string errOut = "EXITING because of fatal PARAMETERS error: --soloCellFilter "
                                 "CellRanger2.2 requires exactly 3 numerical parameters";
            errOut += "\nSOLUTION: re-run with --soloCellFilter CellRanger2.2 <nExpectedCells> "
                      "<maxPercentile> <maxMinRatio>\n";
            exitWithError(errOut, std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER,
                          *context);
        };

        context->inOut->logMain
            << "ParametersSolo: --soloCellFilter CellRanger2.2 filtering parameters: " << pars1
            << std::endl;
        std::istringstream parsStream(pars1);
        parsStream >> cellFilter.knee.nExpectedCells >> cellFilter.knee.maxPercentile >>
            cellFilter.knee.maxMinRatio;

    } else if (cellFilter.type[0] == "EmptyDrops_CR") {
        if (cellFilter.type.size() == 1) {
            context->inOut->logMain << "ParametersSolo: using hardcoded filtering parameters for "
                                       "--soloCellFilter EmptyDrops_CR\n";
            pars1 = "3000 0.99 10 45000 90000 500 0.01 20000 0.01 10000";
        } else if (cellFilter.type.size() < 11) {
            std::string errOut = "EXITING because of fatal PARAMETERS error: --soloCellFilter "
                                 "EmptyDrops_CR requires exactly 10 numerical parameters";
            errOut += "\nSOLUTION: re-run with --soloCellFilter EmptyDrops_CR ";
            errOut += "<nExpectedCells> <maxPercentile> <maxMinRatio> <indMin> <indMax> <umiMin> "
                      "<umiMinFracMedian> <candMaxN> <FDR> <simN>\n";
            exitWithError(errOut, std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER,
                          *context);
        };

        context->inOut->logMain
            << "ParametersSolo: --soloCellFilter EmptyDrops_CR filtering parameters: " << pars1
            << std::endl;
        std::istringstream parsStream(pars1);
        parsStream >> cellFilter.knee.nExpectedCells >> cellFilter.knee.maxPercentile >>
            cellFilter.knee.maxMinRatio;
        parsStream >> cellFilter.eDcr.indMin >> cellFilter.eDcr.indMax >> cellFilter.eDcr.umiMin >>
            cellFilter.eDcr.umiMinFracMedian;
        parsStream >> cellFilter.eDcr.candMaxN >> cellFilter.eDcr.FDR >> cellFilter.eDcr.simN;

    } else if (cellFilter.type[0] == "TopCells") {
        if (cellFilter.type.size() < 2) {
            std::string errOut = "EXITING because of fatal PARAMETERS error: number of cells not "
                                 "specified for --soloCellFilter TopCells";
            errOut += "\nSOLUTION: --soloCellFilter TopCells <NumberOfCells>\n";
            exitWithError(errOut, std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER,
                          *context);
        };
        cellFilter.topCells = stoi(cellFilter.type[1]);
    } else if (cellFilter.type[0] == "None") {
        // nothing to do
    } else {
        exitWithError("EXITING because of fatal PARAMETERS error: unrecognized option in "
                      "--soloCellFilter=" +
                          cellFilter.type[0] +
                          "\nSOLUTION: use allowed options: CellRanger2.2 or None\n",
                      std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER, *context);
    };
};

} // namespace nustar::star2711b::solo::detail
