#include "sequence.hpp"
#include "../configuration/context.hpp"
#include "../configuration/options.hpp"
#include "../support/algorithms.hpp"
#include "../support/error.hpp"
#include "../support/streams.hpp"
#include <stdlib.h>

namespace nustar::star2711b::solo::detail {
void Options::complexWLstrings() {

    cbWLstr.resize(cbWLsize);

    for (auto& cb : cbV) { // initialize
        cb.iCB = 0;
        cb.iLen = cb.minLen;
    };

    for (std::uint32_t ii = 0; ii < cbWLsize; ii++) { // cycle over full WL
        for (UnsignedWord ii = 0; ii < cbV.size();
             ii++) { // check for overflow and re-calculate all indexes
            BarcodeLayout& cb = cbV[ii];
            if (cb.iCB == cb.wl[cb.iLen].size()) { // advance length
                cb.iLen++;
                cb.iCB = 0; // reset iCB
            };
            if (cb.iLen == cb.wl.size()) { // advance barcode
                cbV[ii + 1].iCB++;
                cb.iLen = cb.minLen; // reset length
            };
        };

        for (auto& cb : cbV)
            cbWLstr[ii] += convertNuclInt64toString(cb.wl[cb.iLen][cb.iCB], cb.iLen) + "_";
        cbWLstr[ii].pop_back();

        cbV[0].iCB++; // shift by one for the next CB
    };
};

void Options::loadWhitelist() {
    if (type == SoloTypes::CB_UMI_Simple || type == SoloTypes::CB_samTagOut) { // simple whitelist
        if (soloCBwhitelist.size() > 1) {
            std::ostringstream errOut;
            errOut
                << "EXITING because of FATAL ERROR in INPUT parameters: --soloCBwhitelist contains "
                   "more than one file which is not allowed with --soloType CB_UMI_Simple \n";
            errOut << "SOLUTION: in --soloCBwhitelist specify only one whitelist file \n";
            exitWithError(errOut.str(), std::cerr, context->inOut->logMain, EXIT_CODE_INPUT_FILES,
                          *context);
        } else if (soloCBwhitelist[0] == "-") {
            std::ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR in INPUT parameters: --soloCBwhitelist is "
                      "not defined\n";
            errOut << "SOLUTION: in --soloCBwhitelist specify path to and name of the whitelist "
                      "file, or None for CB demultiplexing without whitelist \n";
            exitWithError(errOut.str(), std::cerr, context->inOut->logMain, EXIT_CODE_INPUT_FILES,
                          *context);
        } else if (soloCBwhitelist[0] == "None") {
            cbWLyes = false;
        } else {
            cbWLyes = true;
            std::ifstream& cbWlStream =
                ifstrOpen(soloCBwhitelist[0], sourceLocation(),
                          "SOLUTION: check the path and permissions of the CB whitelist file: " +
                              soloCBwhitelist[0],
                          *context);
            std::string seq1;
            while (cbWlStream >> seq1) {
                if (seq1.size() != cbL) {
                    std::ostringstream errOut;
                    errOut << "EXITING because of FATAL ERROR in input CB whitelist file: "
                           << soloCBwhitelist[0] << " the total length of barcode sequence is "
                           << seq1.size() << " not equal to expected " << bL << "\n";
                    errOut << "SOLUTION: make sure that the barcode read is the second in "
                              "--readFilesIn and check that is has the correct formatting\n";
                    exitWithError(errOut.str(), std::cerr, context->inOut->logMain,
                                  EXIT_CODE_INPUT_FILES, *context);
                };
                UnsignedWord cb1;
                if (convertNuclStrToInt64(seq1, cb1)) { // convert to 2-bit format
                    cbWL.push_back(cb1);
                } else {
                    context->inOut->logMain
                        << "WARNING: CB whitelist sequence contains non-ACGT base and is ignored: "
                        << seq1 << std::endl;
                };
            };
            if (cbWL.size() == 0) { // empty whitelist
                exitWithError("EXITING because of FATAL ERROR: CB whitelist file " +
                                  soloCBwhitelist[0] +
                                  " is empty. \nSOLUTION: provide non-empty whitelist.\n",
                              std::cerr, context->inOut->logMain, EXIT_CODE_INPUT_FILES, *context);
            };
        };

        std::sort(cbWL.begin(), cbWL.end());              // sort
        auto un1 = std::unique(cbWL.begin(), cbWL.end()); // collapse identical
        cbWL.resize(std::distance(cbWL.begin(), un1));
        cbWLsize = cbWL.size();
        context->inOut->logMain << "Number of CBs in the whitelist = " << cbWLsize << std::endl;

        cbWLstr.resize(cbWLsize);
        for (UnsignedWord ii = 0; ii < cbWLsize; ii++)
            cbWLstr[ii] = convertNuclInt64toString(cbWL[ii], cbL);
    } else if (type == SoloTypes::SmartSeq) {
        cbWLstr = context->outSAMattrRG;
        cbWLsize = cbWLstr.size();
        cbWLyes = true;
    } else if (type == SoloTypes::CB_UMI_Complex) { // complex barcodes: multiple whitelist (one for
                                                    // each CB), varying CB length
        cbWLyes = true; // for complex barcodes, no-whitelist option is not allowed for now

        adapterYes = false;
        if (adapterSeq != "-")
            adapterYes = true;

        if (cbPositionStr.size() != soloCBwhitelist.size()) {
            std::ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETER error: number of barcodes in "
                      "--soloCBposition : "
                   << cbPositionStr.size()
                   << " is not equal to the number of WhiteLists in --soloCBwhitelist : "
                   << soloCBwhitelist.size() << "\n";
            errOut << "SOLUTION: make sure that the number of CB whitelists and CB positions are "
                      "the same\n";
            exitWithError(errOut.str(), std::cerr, context->inOut->logMain, EXIT_CODE_INPUT_FILES,
                          *context);
        };
        cbV.resize(cbPositionStr.size());
        for (std::uint32_t ii = 0; ii < cbPositionStr.size(); ii++) {
            cbV[ii].extractPositionsFromString(cbPositionStr[ii]);
        };

        umiV.extractPositionsFromString(umiPositionStr);
        umiL = 0; // this will be defined when the first barcode is processed

        umiV.adapterLength = adapterSeq.size(); // one adapter for all
        cbWLsize = 1;
        for (std::uint32_t icb = 0; icb < cbV.size(); icb++) { // cycle over WL files
            cbV[icb].adapterLength = adapterSeq.size();        // one adapter for all

            std::ifstream& cbWlStream =
                ifstrOpen(soloCBwhitelist[icb], sourceLocation(),
                          "SOLUTION: check the path and permissions of the CB whitelist file: " +
                              soloCBwhitelist[icb],
                          *context);

            std::string seq1;
            while (cbWlStream >> seq1) { // cycle over one WL file
                UnsignedWord cb1;
                if (!convertNuclStrToInt64(seq1, cb1)) { // convert to 2-bit format
                    context->inOut->logMain
                        << "WARNING: CB whitelist sequence contains non-ACGT base and is ignored: "
                        << seq1 << std::endl;
                    continue;
                };

                std::uint32_t len1 = seq1.size();
                if (len1 >= cbV[icb].wl.size())
                    cbV[icb].wl.resize(len1 + 1); // add new possible lengths to this CB
                cbV[icb].wl.at(len1).push_back(cb1);
            };

            cbV[icb].sortWhiteList(this);
            cbV[icb].wlFactor = cbWLsize;
            cbWLsize *= cbV[icb].totalSize;
        };

        complexWLstrings();
    };

    time_t rawTime;
    time(&rawTime);
    context->inOut->logMain
        << timeMonthDayTime(rawTime)
        << " ... Finished reading, sorting and deduplicating CB whitelist sequences." << std::endl;
}

} // namespace nustar::star2711b::solo::detail
