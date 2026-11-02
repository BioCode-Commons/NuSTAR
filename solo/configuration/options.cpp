#include "options.hpp"
#include "../barcode/sequence.hpp"
#include "context.hpp"
#include "../support/algorithms.hpp"
#include "../support/error.hpp"
#include "../support/streams.hpp"
#include <stdlib.h>

namespace nustar::star2711b::solo::detail {

void Options::initialize(Context* pPin) {
    context = pPin;

    cellFiltering();
    if (context->runMode ==
        "soloCellFiltering") { // only filtering happens, do not need any other parameters
        yes = true;
        umiDedup.typesIn = {"NoDedup"}; // this does not affect the results - the dedup had been
                                        // done when raw matrix was generated
        umiDedup.initialize(this);
        return;
    };

    // constants - may turn into parameters in the future
    redistrReadsNfiles = 3 * context->runThreadN;
    barcodeStart = barcodeEnd = 0; // this means that the entire barcodeRead is considered barcode.
                                   // Will change it for simple barcodes.
    yes = true;
    if (typeStr == "None" || typeStr == "SmartSeq") { // solo SAM attributes not allowed

        if (context->outSAMattrPresent.CR || context->outSAMattrPresent.CY ||
            context->outSAMattrPresent.UR || context->outSAMattrPresent.UY ||
            context->outSAMattrPresent.CB || context->outSAMattrPresent.UB ||
            context->outSAMattrPresent.sS || context->outSAMattrPresent.sQ ||
            context->outSAMattrPresent.sM || context->outSAMattrPresent.sF) {
            std::ostringstream errOut;
            errOut << "EXITING because of FATAL INPUT ERROR: --outSAMattributes contains "
                      "CR/CY/UR/UY/CB/UB tags which are not allowed for --soloType "
                   << typeStr << '\n';
            errOut << "SOLUTION: re-run STAR without these attribures\n";
            exitWithError(errOut.str(), std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER,
                          *context);
        };
    };

    if (typeStr == "None") {
        type = SoloTypes::None;
        yes = false;
        samAttrYes = false;
        return;
    } else if (typeStr == "CB_UMI_Simple" || typeStr == "Droplet") {
        type = SoloTypes::CB_UMI_Simple;
        if (umiL > 16) {
            std::ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: UMI length is too long: "
                      "--soloUMIlen="
                   << umiL << "\n";
            errOut << "SOLUTION: UMI length cannot be longer than 16";
            exitWithError(errOut.str(), std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER,
                          *context);
        };
        if (cbL > 31) {
            std::ostringstream errOut;
            errOut
                << "EXITING because of fatal PARAMETERS error: CB length is too long: --soloCBlen="
                << cbL << "\n";
            errOut << "SOLUTION: CB length cannot be longer than 31";
            exitWithError(errOut.str(), std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER,
                          *context);
        };

        cbumiL = cbL + umiL;
        if (bL == 1)
            bL = cbumiL;

        barcodeStart = std::min(cbS, umiS) - 1;
        barcodeEnd = std::max(cbS + cbL, umiS + umiL) - 2;

    } else if (typeStr == "CB_UMI_Complex") {
        type = SoloTypes::CB_UMI_Complex;
        bL = 0; // no fixed length for barcode sequence
        cbumiL = 0;

    } else if (typeStr == "CB_samTagOut") {
        type = SoloTypes::CB_samTagOut;
        cbumiL = cbL + umiL;
        if (bL == 1)
            bL = cbumiL;

        barcodeStart = std::min(cbS, umiS) - 1;
        barcodeEnd = std::max(cbS + cbL, umiS + umiL) - 2;

    } else if (typeStr == "SmartSeq") {
        type = SoloTypes::SmartSeq;

    } else {
        std::ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETERS error: unrecognized option in --soloType="
               << typeStr << "\n";
        errOut << "SOLUTION: use allowed option: None OR CB_UMI_Simple OR CB_UMI_Complex\n";
        errOut << "Obsolete option Droplet should be replaced with CB_UMI_Simple";
        exitWithError(errOut.str(), std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER,
                      *context);
    };
    if (CBtype.typeString == "Sequence") {
        CBtype.type = 1;
    } else if (CBtype.typeString == "String") {
        if (soloCBwhitelist[0] != "None") {
            std::ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: --soloCBtype String cannot be "
                      "used with passlist\n";
            errOut << "SOLUTION: use --soloCBwhitelist None\n";
            exitWithError(errOut.str(), std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER,
                          *context);
        };
        CBtype.type = 2;
        CBtype.strMtx = std::make_unique<std::mutex>();
        CBtype.strMap.reserve(100000);
    } else {
        std::ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETERS error: unrecognized option in --soloCBtype="
               << CBtype.typeString << "\n";
        errOut << "SOLUTION: use allowed option: Sequence OR String\n";
        exitWithError(errOut.str(), std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER,
                      *context);
    };
    barcodeRead = -1;
    barcodeReadSeparate = false;

    if (context->readFilesTypeN != 10) {   // input from FASTQ
        if (type == SoloTypes::SmartSeq) { // no barcode read
            // TODO: a lot of parameters should not be defined for SmartSeq option - check it here
        } else { // all other types require barcode read
            if (barcodeReadIn ==
                0) { // barcode read is separate - needs to be the last read in --readFilesIn
                if (context->readNends < 2) {
                    exitWithError(
                        "EXITING because of fatal PARAMETERS error: --soloType (except SmartSeq) "
                        "with --soloBarcodeMate 0 (default) require 2 reads or 3 reads, where the "
                        "last read is the barcode read.\n"
                        "SOLUTION: if barcode is in a separate mate, specify it as the last file "
                        "in --readFilesIn. If barcode sequence is a part of one of the mates, "
                        "specify that mate with --soloBarcodeMate 1 (or 2 or 3)",
                        std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER, *context);
                };
                context->readNmates = context->readNends - 1; // true mates, excluding barcode read
                barcodeRead = context->readNends - 1; // the barcode read is always the last one
            } else if (barcodeReadIn > context->readNends) {
                exitWithError(
                    "EXITING because of fatal PARAMETERS error: --soloBarcodeMate " +
                        std::to_string(barcodeReadIn) + "is larger than number of mates " +
                        std::to_string(context->readNends) +
                        "\nSOLUTION: specify --soloBarcodeMate <= than the number of mates.",
                    std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER, *context);
            } else { // barcode sequence is in one of the mates
                if (type != SoloTypes::CB_UMI_Simple) {
                    exitWithError("EXITING because of fatal PARAMETERS error: --soloBarcodeMate " +
                                      std::to_string(barcodeReadIn) +
                                      ">0 for is not allowed for --soloType " + typeStr +
                                      "\nSOLUTION: specify --soloBarcodeMate 0   or   --soloType "
                                      "CB_UMI_Simple",
                                  std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER,
                                  *context);
                };

                barcodeRead = barcodeReadIn - 1;
                barcodeReadSeparate = true;
                bL = 0; //>0 only for separate barcode read

                if (context->pClip.in[0].N[barcodeRead] == 0 &&
                    context->pClip.in[1].N[barcodeRead] ==
                        0) { // clipping not specified for the mate with barcodes
                    exitWithError(
                        "EXITING because of fatal PARAMETERS error: --soloBarcodeMate " +
                            std::to_string(barcodeReadIn) +
                            " specifies that barcode sequence is a part of the mate " +
                            std::to_string(barcodeReadIn) +
                            ", which requires clipping the barcode off this mate."
                            "\nSOLUTION: clip the barcode sequence from 5' or/and 3' with "
                            "--clip5pNbases   or/and --clip3pNbases . The values for mate1 and "
                            "mate2 have to specified, specify 0 for no clipping.",
                        std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER, *context);
                };
            };
        };
    } else if (context->readFilesTypeN == 10) { // input from SAM
        if (typeStr == "SmartSeq") {
            exitWithError(
                "EXITING because of fatal PARAMETERS error: --readFilesType SAM SE/PE cannot be "
                "used with --soloType SmartSeq\n"
                "SOLUTION: for Smart-seq input from BAM files, use --soloType CB_UMI_Simple , "
                "create whitelist of SmartSeq file names, and specify the SAM tag that records "
                "these file names in --soloInputSAMattrBarcodeSeq",
                std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER, *context);
        };

        if (samAtrrBarcodeSeq.at(0) == "-") {
            exitWithError("EXITING because of fatal PARAMETERS error: --readFilesType SAM SE/PE "
                          "requires --soloInputSAMattrBarcodeSeq.\n"
                          "SOLUTION: specify input SAM attributes for barcode sequence in "
                          "--soloInputSAMattrBarcodeSeq, and (optionally) quality with "
                          "--soloInputSAMattrBarcodeQual",
                          std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER, *context);
        };

        if (samAtrrBarcodeQual.at(0) == "-") {
            warningMessage(" since --readFilesType SAM SE/PE --soloInputSAMattrBarcodeQual - : "
                           "qualities for barcode read will be replaced with 'H'",
                           context->inOut->logMain, std::cerr, *context);
            samAtrrBarcodeQual.clear();
        };

        for (auto& tag : samAtrrBarcodeSeq) {
            if (tag.size() != 2) {
                exitWithError(
                    "EXITING because of fatal PARAMETERS error: --soloInputSAMattrBarcodeSeq "
                    "attributes have to be two-letter strings.\n"
                    "SOLUTION: specify correct two-letter strings in --soloInputSAMattrBarcodeSeq",
                    std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER, *context);
            };
            tag.insert(0, 1, '\t');
        };

        for (auto& tag : samAtrrBarcodeQual) {
            if (tag.size() != 2) {
                exitWithError(
                    "EXITING because of fatal PARAMETERS error: --soloInputSAMattrBarcodeQual "
                    "attributes have to be two-letter strings.\n"
                    "SOLUTION: specify correct two-letter strings in --soloInputSAMattrBarcodeQual",
                    std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER, *context);
            };
            tag.insert(0, 1, '\t');
        };
    };
    if (strandStr == "Unstranded") {
        strand = -1;
    } else if (strandStr == "Forward") {
        strand = 0;
    } else if (strandStr == "Reverse") {
        strand = 1;
    } else {
        std::ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETERS error: unrecognized option in --soloStrand="
               << strandStr << "\n";
        errOut << "SOLUTION: use allowed option: Unstranded OR Forward OR Reverse";
        exitWithError(errOut.str(), std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER,
                      *context);
    };
    featureInd.fill(-1);
    featureYes.fill(false);
    featureFirst = -1;
    for (auto& fin : featureIn) {
        for (std::uint32_t ii = 0; ii < FeatureTypes::Names.size(); ii++) {
            if (fin == FeatureTypes::Names[ii]) {
                featureYes[ii] = true;
                features.push_back(ii);
                if (featureFirst == -1)
                    featureFirst = ii; // recored the first feature on the list
                break;
            };
        };
    };
    nFeatures = features.size();
    std::sort(features.begin(), features.end());
    for (std::uint32_t ii = 0; ii < nFeatures; ii++)
        featureInd[features[ii]] = ii;

    if (nFeatures != featureIn.size()) {
        std::ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETERS error: unrecognized option(s) in "
                  "--soloFeatures ";
        for (auto& fname : featureIn)
            errOut << fname << " ";
        errOut << "\nSOLUTION: use allowed option: ";
        errOut << FeatureTypes::Names[0];
        for (auto& fname : FeatureTypes::Names)
            errOut << "   OR   " << fname;
        exitWithError(errOut.str(), std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER,
                      *context);
    };

    if (featureYes[FeatureTypes::Velocyto] && type == SoloTypes::SmartSeq) {
        std::string errOut = "EXITING because of fatal PARAMETERS error: --soloFeatures Velocyto "
                             "is presently not compatible with --soloType SmartSeq .\n";
        errOut += "SOLUTION: re-run without --soloFeatures Velocyto .";
        exitWithError(errOut, std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER, *context);
    };

    if (featureYes[FeatureTypes::Gene]) {
        context->quant.gene.yes = true;
        context->quant.yes = true;
    };
    if (featureYes[FeatureTypes::GeneFull]) {
        context->quant.geneFull.yes = true;
        context->quant.yes = true;
        if (!featureYes[FeatureTypes::Gene])
            context->quant.gene.yes =
                false; // if GeneFull is requested, but Gene is not, turn it off - it could have
                       // been turned on because of GX/GN attributes
    };
    if (featureYes[FeatureTypes::GeneFull_Ex50pAS]) {
        context->quant.geneFull_Ex50pAS.yes = true;
        context->quant.yes = true;
        if (!featureYes[FeatureTypes::Gene])
            context->quant.gene.yes =
                false; // if GeneFull is requested, but Gene is not, turn it off - it could have
                       // been turned on because of GX/GN attributes
    };
    if (featureYes[FeatureTypes::GeneFull_ExonOverIntron]) {
        context->quant.geneFull_ExonOverIntron.yes = true;
        context->quant.gene.yes = true; // needed to prioritize exons over introns
        context->quant.yes = true;
    };

    // initialize CB match to WL types
    init_CBmatchWL();
    umiDedup.initialize(this);

    // make output directory if needed
    if (outFileNames[0].find_last_of("/") < outFileNames[0].size()) { // need to create dir
        std::string dir1 = context->outFileNamePrefix +
                           outFileNames[0].substr(0, outFileNames[0].find_last_of("/"));
        if (mkdir(dir1.c_str(), context->runDirPerm) != 0 && errno != EEXIST) {
            std::ostringstream errOut;
            errOut << "EXITING because of fatal OUTPUT FILE error: could not create Solo output "
                      "directory"
                   << dir1 << "\n";
            errOut << "SOLUTION: check the path and permisssions";
            exitWithError(errOut.str(), std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER,
                          *context);
        };
    };

    QSbase = 33; // TODO make these user-definable
    QSmax = 33;
    cbMinP = 0.975;

    umiMaskLow = (std::uint32_t)((((UnsignedWord)1) << umiL) - 1);
    umiMaskHigh = ~umiMaskLow;

    loadWhitelist();
    samAttrYes = false;
    if ((context->outSAMattrPresent.CB || context->outSAMattrPresent.UB) &&
        type != SoloTypes::CB_samTagOut) {
        samAttrYes = true;
        if (!context->outBAMcoord) {
            std::ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: CB and/or UB attributes in "
                      "--outSAMattributes can only be output in the sorted BAM file.\n";
            errOut << "SOLUTION: re-run STAR with --outSAMtype BAM SortedByCoordinate ...\n";
            exitWithError(errOut.str(), std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER,
                          *context);
        };
    } else if (context->outSAMattrPresent.UB && type == SoloTypes::CB_samTagOut) {
        exitWithError("EXITING because of fatal PARAMETERS error: UB attribute (corrected UMI) in "
                      "--outSAMattributes cannot be used with --soloType CB_samTagOut \n"
                      "SOLUTION: instead, use UR (uncorrected UMI) in --outSAMattributes\n",
                      std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER, *context);
    };
    readInfoYes.fill(false);
    if (featureYes[FeatureTypes::VelocytoSimple] ||
        featureYes[FeatureTypes::Velocyto]) { // turn readInfo on for Gene needed by VelocytoSimple
        readInfoYes[FeatureTypes::Gene] = true;
    };
    samAttrFeature = featureFirst;
    if (samAttrYes) { // pSolo.samAttrFeature=0 by default, so need to check samAttrYes
        if (featureFirst == FeatureTypes::Gene || featureFirst == FeatureTypes::GeneFull ||
            featureFirst == FeatureTypes::GeneFull_Ex50pAS ||
            featureFirst == FeatureTypes::GeneFull_ExonOverIntron) {
            // all good
        } else {
            std::ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: CB and/or UB attributes in "
                      "--outSAMattributes require --soloFeatures Gene OR/AND GeneFull OR/AND "
                      "GeneFull_Ex50pAS.\n";
            errOut << "SOLUTION: re-run STAR adding Gene AND/OR GeneFull OR/AND GeneFull_Ex50pAS "
                      "OR/AND GeneFull_ExonOverIntron to --soloFeatures\n";
            exitWithError(errOut.str(), std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER,
                          *context);
        };
        readInfoYes[samAttrFeature] = true;
    };
    readIndexYes = readInfoYes;
    readStats.yes = false;
    readStatsYes.fill(false);
    if (readStats.type == "Standard") {
        readStats.yes = true;
        readStatsYes.fill(true);
        readStatsYes[FeatureTypes::VelocytoSimple] =
            false; // this could be allowed, but it will have the same info as Gene
        readStatsYes[FeatureTypes::Velocyto] = false;
        readStatsYes[FeatureTypes::SJ] =
            false; // output for SJ requires careful consideration for SoloReadFeature_record.cpp,
                   // for output of reads that do not have splice junction
        for (std::uint32_t ff = 0; ff < readIndexYes.size(); ff++) { // merge with previous values
            readIndexYes[ff] |= readStatsYes[ff];
        };
    } else if (readStats.type != "None") {
        exitWithError("EXITING because of fatal PARAMETERS error: unrecognized option in "
                      "--soloCellReadStats" +
                          readStats.type + "\nSOLUTION: use allowed options: None OR Standard \n",
                      std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER, *context);
    };
    if (umiFiltering.type[0] == "MultiGeneUMI") {
        umiFiltering.MultiGeneUMI = true;
        umiFiltering.yes = true;
    } else if (umiFiltering.type[0] == "MultiGeneUMI_All") {
        umiFiltering.MultiGeneUMI_All = true;
        umiFiltering.yes = true;
    } else if (umiFiltering.type[0] == "MultiGeneUMI_CR") {
        umiFiltering.MultiGeneUMI_CR = true;
        if (umiDedup.typesIn.size() > 1 || umiDedup.typesIn.at(0) != "1MM_CR")
            exitWithError("EXITING because of fatal PARAMETERS error: --soloUMIfiltering "
                          "MultiGeneUMI_CR only works with --soloUMIdedup 1MM_CR"
                          "\nSOLUTION: rerun with --soloUMIfiltering MultiGeneUMI_CR "
                          "--soloUMIdedup 1MM_CR \n",
                          std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER, *context);
    } else if (umiFiltering.type[0] == "-") {
        // nothing to do
    } else {
        exitWithError(
            "EXITING because of fatal PARAMETERS error: unrecognized option in "
            "--soloUMIfiltering=" +
                umiFiltering.type[0] +
                "\nSOLUTION: use allowed options: - or MultiGeneUMI or MultiGeneUMI_CR \n",
            std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER, *context);
    };
    multiMap.initialize(this);
    if (multiMap.yes.multi) {
        if (type == SoloTypes::CB_samTagOut || type == SoloTypes::SmartSeq) {
            exitWithError("EXITING because of fatal PARAMETERS error: multimapping options do not "
                          "work for --soloType " +
                              typeStr +
                              "\nSOLUTION: use default option --soloMultiMappers Unique\n",
                          std::cerr, context->inOut->logMain, EXIT_CODE_PARAMETER, *context);
        };
        readIndexYes[FeatureTypes::Gene] = true;
        readIndexYes[FeatureTypes::GeneFull] = true;
        readIndexYes[FeatureTypes::GeneFull_Ex50pAS] = true;
        readIndexYes[FeatureTypes::GeneFull_ExonOverIntron] = true;
    };
};

} // namespace nustar::star2711b::solo::detail
