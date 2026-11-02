#pragma once

#include "../annotation/feature_types.hpp"
#include "../barcode/layout.hpp"
#include "../support/types.hpp"
#include <array>
#include <memory>
#include <mutex>
#include <unordered_map>

namespace nustar::star2711b::solo::detail {

class Context;
class Options;

class UmiDedupOptions {
public:
    const static std::uint32_t tN = 6;
    std::array<std::string, tN> typeNames{
        {"NoDedup", "Exact", "1MM_All", "1MM_Directional", "1MM_CR", "1MM_Directional_UMItools"}};
    enum typeI : std::int32_t {
        NoDedup = 0,
        Exact = 1,
        All = 2,
        Directional = 3,
        CR = 4,
        Directional_UMItools = 5
    };

    struct {
        uint32_t N;
        std::array<bool, tN> B;
        bool &NoDedup = B[0], &Exact = B[1], &All = B[2], &Directional = B[3], &CR = B[4],
             &Directional_UMItools = B[5];
    } yes;

    struct {
        // uint32_t N;
        std::array<uint32_t, tN> I;
        uint32_t &NoDedup = I[0], &Exact = I[1], &All = I[2], &Directional = I[3], &CR = I[4],
                 &Directional_UMItools = I[5];
        uint32_t main; // index for SAM/stats/filtering output
    } countInd;        // index in the countCellGennUMI

    std::vector<std::string> typesIn; // UMIdedup types from user options
    std::vector<std::int32_t> types;  // the above converted to typeI numbers
    std::int32_t
        typeMain; // the type to be used in SAM/stats/filtering output - for now just types[0]

    void initialize(Options* pS);
    //    int it;
};

class MultimapperOptions {
public:
    const static std::uint32_t tN = 5;
    std::array<std::string, tN> typeNames{{"Unique", "Uniform", "Rescue", "PropUnique", "EM"}};
    enum typeI : std::int32_t { Unique = 0, Uniform = 1, Rescue = 2, PropUnique = 3, EM = 4 };

    struct {
        bool multi; // if multimappers are requested
        uint32_t N;
        std::array<bool, tN> B;
        bool &Unique = B[0], &Uniform = B[1], &Rescue = B[2], &PropUnique = B[3], &EM = B[4];
    } yes;

    struct {
        // uint32_t N;
        std::array<uint32_t, tN> I;
        uint32_t &Unique = I[0], &Uniform = I[1], &Rescue = I[2], &PropUnique = I[3], &EM = I[4];
        uint32_t main; // index for SAM/stats/filtering output
    } countInd;        // index in the countCellGennUMI

    std::vector<std::string> typesIn; // UMIdedup types from user options
    std::vector<std::int32_t> types;  // the above converted to typeI numbers
    std::int32_t
        typeMain; // the type to be used in SAM/stats/filtering output - for now just types[0]

    void initialize(Options* pS);
};

class Options {
public:
    Context* context;
    bool yes;

    // chemistry, library etc
    std::string typeStr;
    enum SoloTypes : std::int32_t {
        None = 0,
        CB_UMI_Simple = 1,
        CB_UMI_Complex = 2,
        CB_samTagOut = 3,
        SmartSeq = 4
    };
    SoloTypes type;
    std::string strandStr;
    std::int32_t strand;

    std::uint32_t barcodeRead, barcodeReadIn; // which read is the barcode read = 0,1,2?
    std::uint32_t barcodeStart, barcodeEnd;   // start/end of barcode sequence on barcodeRead
    bool barcodeReadSeparate;

    // simple barcodes
    std::uint32_t cbS, cbL;              // cell barcode start,length
    std::uint32_t umiS, umiL;            // umi start,length
    std::once_flag complexUmiLengthOnce; // first extracted complex UMI, shared by workers
    std::uint32_t bL, cbumiL; // total barcode sequene length, CB+UMI length. Former does may not be
                              // equal to the latter

    std::vector<std::string> cbPositionStr;
    std::string umiPositionStr;

    // complex barcodes
    std::vector<BarcodeLayout> cbV;
    BarcodeLayout umiV;                  // single UMI
    bool adapterYes;                     // anchor?
    std::string adapterSeq;              // anchor sequence
    std::uint32_t adapterMismatchesNmax; // max number of mismatches in the anchor

    // input from SAM files
    std::vector<std::string> samAtrrBarcodeSeq, samAtrrBarcodeQual;

    // whitelist - general
    UnsignedWord cbWLsize;
    bool cbWLyes;
    std::vector<std::string> soloCBwhitelist;
    std::vector<UnsignedWord> cbWL;
    std::vector<std::string> cbWLstr;

    MultimapperOptions multiMap;

    // features
    std::vector<std::string> featureIn; // string of requested features
    std::vector<std::uint32_t> features;
    std::uint32_t nFeatures; //=features.size(), number of requested features

    std::int32_t featureFirst;                    // which feature is the first on the list
    std::array<bool, FeatureTypes::N> featureYes; // which features are requested
    std::array<bool, FeatureTypes::N>
        readInfoYes; // which features will need readInfo (for now only Gene and GeneFull)
    std::array<bool, FeatureTypes::N>
        readIndexYes; // which features will need recording of readIndex (for now only Gene and
                      // GeneFull, for multimappers)
    std::array<bool, FeatureTypes::N>
        readStatsYes; // which features will need output of read statistics
    std::array<std::int32_t, FeatureTypes::N>
        featureInd; // index of each feature - skips unrequested features

    // filtering
    char QSbase, QSmax; // quality score base and cutoff

#ifdef MATCH_CellRanger
    double cbMinP; // for CBs with non-exact matching to WL, min posterior probability
#else
    float cbMinP; // for CBs with non-exact matching to WL, min posterior probability
#endif

    // cell filtering
    struct {
        std::vector<std::string> type;
        std::uint32_t topCells;

        struct {
            double nExpectedCells;
            double maxPercentile;
            double maxMinRatio;
        } knee;

        struct {
            std::uint32_t indMin, indMax; // min/max cell index, sorted by UMI counts,for empty
                                          // cells
            std::uint32_t umiMin;
            double umiMinFracMedian;
            std::uint32_t candMaxN;
            double FDR;
            std::uint32_t simN;
        } eDcr; // EmptyDrops-CellRanger

    } cellFilter;

    // CBtype
    struct {
        std::string typeString;
        std::int32_t type = 0;
        std::unordered_map<std::string, std::uint32_t> strMap;
        std::unique_ptr<std::mutex> strMtx;
    } CBtype;

    // CB match
    struct {
        std::string type;
        bool mm1;             // 1 mismatch allowed
        bool mm1_multi;       // 1 mismatch, multiple matches to WL allowed
        bool oneExact;        // CBs require at least one exact match
        bool mm1_multi_pc;    // use psedocounts while calculating probabilities of multi-matches
        bool mm1_multi_Nbase; // allow multimatching to WL for CBs with N-bases
        bool EditDist_2;      // allow EditDistance <=3
    } CBmatchWL;

    // UMIdedup
    UmiDedupOptions umiDedup;

    // multi-gene umi
    struct {
        std::vector<std::string> type;
        bool MultiGeneUMI = false;
        bool MultiGeneUMI_All = false;
        bool yes = false; // true for non-CR
        bool MultiGeneUMI_CR = false;
    } umiFiltering;

    // clusters
    std::string clusterCBfile;

    // output
    std::vector<std::string> outFileNames;
    struct {
        std::string featuresGeneField3;
    } outFormat;

    bool samAttrYes;             // post-processed SAM attributes: error-corrected CB and UMI
    std::int32_t samAttrFeature; // which feature to use for error correction

    // processing
    std::uint32_t redistrReadsNfiles; // numer of files to resditribute reads into

    // constants
    std::uint32_t umiMaskLow, umiMaskHigh; // low/high half bit-mask or UMIs

    // readStats
    struct {
        std::string type; // input parameter
        bool yes = false;
    } readStats;

    void initialize(Context* pPin);
    void umiSwapHalves(std::uint32_t& umi);
    void complexWLstrings();
    void loadWhitelist();
    void cellFiltering();

    void init_CBmatchWL();
};

} // namespace nustar::star2711b::solo::detail
