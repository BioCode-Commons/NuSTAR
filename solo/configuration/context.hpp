#pragma once

#include "options.hpp"
#include "parameter.hpp"
#include "../support/log_streams.hpp"
#include "../support/time.hpp"
#include "../support/types.hpp"
#include <memory>
#include <mutex>
#include <unordered_set>

namespace nustar::star2711b::solo::detail {

// Shared run configuration and stream ownership. Mapping/index settings live
// in NuSTAR's runtime; counting receives only the fields it consumes.
class Context {
public:
    Context();
    ~Context() = default;
    Context(const Context&) = delete;
    Context& operator=(const Context&) = delete;
    std::vector<std::unique_ptr<std::ios>> ownedStreams;
    std::mutex streamMutex, errorMutex;
    struct {
        UnsignedWord readN = 0, mappedReadsU = 0, mappedReadsM = 0;
    } countingStats;
    std::vector<std::unique_ptr<ParameterBase>> parArray;
    std::vector<std::string> parameterInputName;
    std::string runMode;
    std::vector<std::string> runModeIn;
    int runThreadN = 1;
    mode_t runDirPerm = 0700;
    struct {
        std::int32_t type = 0;
    } runRestart;
    int readFilesTypeN = 2;
    UnsignedWord readNmates = 1, readNends = 2;
    std::uint32_t readQualityScoreBase = 33;
    struct {
        struct {
            std::vector<std::uint32_t> N;
        } in[2];
    } pClip;
    struct {
        std::string gDir, sjdbGTFfile = "-";
        struct {
            bool outQuant = false;
        } transform;
        struct {
            std::unordered_set<UnsignedWord> mito;
        } chrSet;
    } pGe, pGeOut;
    struct {
        std::string outDir;
    } sjdbInsert;
    std::string outFileNamePrefix, outFileTmp;
    std::vector<std::string> outSAMattrRG;
    struct {
        bool CR = false, CY = false, UR = false, UY = false, CB = false, UB = false, sM = false,
             sS = false, sQ = false, sF = false;
    } outSAMattrPresent;
    bool outBAMcoord = false;
    std::vector<std::string> outputAttributes;
    struct {
        bool yes = false;
        struct {
            bool yes = false;
        } trSAM, geneFull, geneFull_Ex50pAS, geneFull_ExonOverIntron, gene;
        struct {
            bool yes = false;
        } geCount;
    } quant;
    Options options;
    std::array<std::vector<UnsignedWord>, 2> sjAll;
    std::unique_ptr<LogStreams> inOut = std::make_unique<LogStreams>();
    void readOptions(const std::vector<std::string>& arguments);
    int scanOneLine(std::string&, int inputLevel, int inputLevelRequested);
    void scanAllLines(std::istream&, int inputLevel, int inputLevelRequested);
};

} // namespace nustar::star2711b::solo::detail
