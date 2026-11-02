#include "context.hpp"
#include "defaults.hpp"
#include "../support/error.hpp"

namespace nustar::star2711b::solo::detail {

constexpr int PAR_NAME_PRINT_WIDTH = 30;
Context::Context() {
    parArray.push_back(
        std::make_unique<VectorParameter<std::string>>(-1, -1, "runMode", &runModeIn));
    parArray.push_back(std::make_unique<ScalarParameter<int>>(-1, -1, "runThreadN", &runThreadN));
    parArray.push_back(std::make_unique<ScalarParameter<std::uint32_t>>(
        -1, -1, "readQualityScoreBase", &readQualityScoreBase));
    parArray.push_back(
        std::make_unique<ScalarParameter<std::string>>(-1, -1, "soloType", &options.typeStr));
    parArray.push_back(
        std::make_unique<ScalarParameter<std::uint32_t>>(-1, -1, "soloCBstart", &options.cbS));
    parArray.push_back(
        std::make_unique<ScalarParameter<std::uint32_t>>(-1, -1, "soloUMIstart", &options.umiS));
    parArray.push_back(
        std::make_unique<ScalarParameter<std::uint32_t>>(-1, -1, "soloCBlen", &options.cbL));
    parArray.push_back(
        std::make_unique<ScalarParameter<std::uint32_t>>(-1, -1, "soloUMIlen", &options.umiL));
    parArray.push_back(std::make_unique<ScalarParameter<std::uint32_t>>(
        -1, -1, "soloBarcodeReadLength", &options.bL));
    parArray.push_back(std::make_unique<ScalarParameter<std::uint32_t>>(-1, -1, "soloBarcodeMate",
                                                                        &options.barcodeReadIn));
    parArray.push_back(std::make_unique<VectorParameter<std::string>>(-1, -1, "soloCBwhitelist",
                                                                      &options.soloCBwhitelist));
    parArray.push_back(
        std::make_unique<ScalarParameter<std::string>>(-1, -1, "soloStrand", &options.strandStr));
    parArray.push_back(std::make_unique<VectorParameter<std::string>>(-1, -1, "soloOutFileNames",
                                                                      &options.outFileNames));
    parArray.push_back(
        std::make_unique<VectorParameter<std::string>>(-1, -1, "soloFeatures", &options.featureIn));
    parArray.push_back(std::make_unique<VectorParameter<std::string>>(-1, -1, "soloUMIdedup",
                                                                      &options.umiDedup.typesIn));
    parArray.push_back(std::make_unique<ScalarParameter<std::string>>(-1, -1, "soloAdapterSequence",
                                                                      &options.adapterSeq));
    parArray.push_back(std::make_unique<ScalarParameter<std::uint32_t>>(
        -1, -1, "soloAdapterMismatchesNmax", &options.adapterMismatchesNmax));
    parArray.push_back(std::make_unique<ScalarParameter<std::string>>(-1, -1, "soloCBmatchWLtype",
                                                                      &options.CBmatchWL.type));
    parArray.push_back(std::make_unique<VectorParameter<std::string>>(-1, -1, "soloCBposition",
                                                                      &options.cbPositionStr));
    parArray.push_back(std::make_unique<ScalarParameter<std::string>>(-1, -1, "soloUMIposition",
                                                                      &options.umiPositionStr));
    parArray.push_back(std::make_unique<VectorParameter<std::string>>(-1, -1, "soloCellFilter",
                                                                      &options.cellFilter.type));
    parArray.push_back(std::make_unique<VectorParameter<std::string>>(-1, -1, "soloUMIfiltering",
                                                                      &options.umiFiltering.type));
    parArray.push_back(std::make_unique<VectorParameter<std::string>>(-1, -1, "soloMultiMappers",
                                                                      &options.multiMap.typesIn));
    parArray.push_back(std::make_unique<ScalarParameter<std::string>>(-1, -1, "soloClusterCBfile",
                                                                      &options.clusterCBfile));
    parArray.push_back(std::make_unique<ScalarParameter<std::string>>(
        -1, -1, "soloOutFormatFeaturesGeneField3", &options.outFormat.featuresGeneField3));
    parArray.push_back(std::make_unique<VectorParameter<std::string>>(
        -1, -1, "soloInputSAMattrBarcodeSeq", &options.samAtrrBarcodeSeq));
    parArray.push_back(std::make_unique<VectorParameter<std::string>>(
        -1, -1, "soloInputSAMattrBarcodeQual", &options.samAtrrBarcodeQual));
    parArray.push_back(std::make_unique<ScalarParameter<std::string>>(-1, -1, "soloCellReadStats",
                                                                      &options.readStats.type));
    parArray.push_back(std::make_unique<ScalarParameter<std::string>>(-1, -1, "soloCBtype",
                                                                      &options.CBtype.typeString));
    parameterInputName = {"Default", "Command-Line-Initial", "Command-Line",
                          "genomeParameters.txt"};
}

void Context::scanAllLines(std::istream& streamIn, int inputLevel, int inputLevelRequested) { // scan
    //     istringstream stringInStream (stringIn);
    std::string lineIn;
    while (std::getline(streamIn, lineIn)) {
        scanOneLine(lineIn, inputLevel, inputLevelRequested);
    };
};

int Context::scanOneLine(std::string& lineIn, int inputLevel,
                         int inputLevelRequested) { // scan one line and load the parameters,
    // 0 if comment, 1 if OK
    if (lineIn == "")
        return 0; // empty line

    std::istringstream lineInStream(lineIn);

    if (inputLevel == 0 && (lineIn.substr(0, 1) == " " || lineIn.substr(0, 1) == "\t"))
        return 0; // for Default input spaces also mark comments, for nice formatting

    std::string parIn("");
    lineInStream >> parIn;
    if (parIn == "" || parIn.substr(0, 2) == "//" || parIn.substr(0, 1) == "#")
        return 0; // this is a comment

    UnsignedWord iPar;
    for (iPar = 0; iPar < parArray.size(); iPar++) {
        if (parIn == parArray[iPar]->nameString) { //
            if (inputLevelRequested < 0 ||
                inputLevelRequested == parArray[iPar]->inputLevelAllowed) {
                break; // will read this parameter values
            } else {
                return 1; // do not read inputs not requested at this level
            };
        };
    };

    std::string parV("");
    lineInStream >> parV;
    if (parV == "") { // parameter value cannot be empty
        std::ostringstream errOut;
        errOut << "EXITING: FATAL INPUT ERROR: empty value for parameter \"" << parIn
               << "\" in input \"" << parameterInputName.at(inputLevel) << "\"\n";
        errOut << "SOLUTION: use non-empty value for this parameter\n" << std::flush;
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    lineInStream.str(lineIn);
    lineInStream.clear();
    lineInStream >> parIn; // get the correct state of stream, past reading parIn

    if (iPar == parArray.size()) { // string is not identified
        std::ostringstream errOut;
        errOut << "EXITING: FATAL INPUT ERROR: unrecognized parameter name \"" << parIn
               << "\" in input \"" << parameterInputName.at(inputLevel) << "\"\n";
        errOut << "SOLUTION: use correct parameter name (check the manual)\n" << std::flush;
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    } else { // found the corresponding parameter
        if (inputLevel == 0 &&
            parArray[iPar]->inputLevel > 0) { // this is one of the initial parameters, it was read
                                              // from Command Line and should not be re-defined
            std::getline(lineInStream, parV);
            inOut->logMain << std::setiosflags(std::ios::left) << std::setw(PAR_NAME_PRINT_WIDTH)
                           << parArray[iPar]->nameString << parV
                           << " ... is RE-DEFINED on Command Line as: " << *(parArray[iPar])
                           << "\n";
        } else if (parArray[iPar]->inputLevelAllowed > 0 &&
                   parArray[iPar]->inputLevelAllowed <
                       inputLevel) { // this is initial parameter and cannot be redefined
            std::ostringstream errOut;
            errOut << "EXITING: FATAL INPUT ERROR: parameter \"" << parIn
                   << "\" cannot be defined at the input level \""
                   << parameterInputName.at(inputLevel) << "\"\n";
            errOut << "SOLUTION: define parameter \"" << parIn << "\" in \""
                   << parameterInputName.at(parArray[iPar]->inputLevelAllowed) << "\"\n"
                   << std::flush;
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        } else if (parArray[iPar]->inputLevel ==
                   inputLevel) { // this parameter was already defined at this input level
            std::ostringstream errOut;
            errOut << "EXITING: FATAL INPUT ERROR: duplicate parameter \"" << parIn
                   << "\" in input \"" << parameterInputName.at(inputLevel) << "\"\n";
            errOut
                << "SOLUTION: keep only one definition of input parameters in each input source\n"
                << std::flush;
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        } else { // read values
            parArray[iPar]->inputValues(lineInStream);
            parArray[iPar]->inputLevel = inputLevel;
            if (inOut->logMain.good()) {
                inOut->logMain << std::setiosflags(std::ios::left)
                               << std::setw(PAR_NAME_PRINT_WIDTH) << parArray[iPar]->nameString
                               << *(parArray[iPar]);
                if (parArray[iPar]->inputLevel > 0)
                    inOut->logMain << "     ~RE-DEFINED";
                inOut->logMain << std::endl;
            };
        };
    };
    return 0;
};

void Context::readOptions(const std::vector<std::string>& arguments) {
    std::istringstream defaults(solo_defaults);
    scanAllLines(defaults, 0, -1);
    std::ostringstream overrides;
    for (const auto& token : arguments) {
        if (token.starts_with("--"))
            overrides << '\n' << token.substr(2);
        else if (token.find_first_of(" \t") != std::string::npos)
            overrides << " \"" << token << '"';
        else
            overrides << ' ' << token;
    }
    std::istringstream options(overrides.str());
    scanAllLines(options, 2, -1);
}

} // namespace nustar::star2711b::solo::detail
