#pragma once

#include "statistics.hpp"
#include "../configuration/context.hpp"
#include "../support/types.hpp"
#include <set>

namespace nustar::star2711b::solo::detail {

class ReadBarcode {
private:
    Context& context_;

public:
    Options& options;

    std::uint32_t homoPolymer[4]; // homopolymer constants
    std::string cbSeq, umiSeq, cbQual, umiQual, bSeq, bQual;
    std::vector<std::string> bStrings; // barcode strings from SAM tags
    std::string cbSeqCorrected;
    UnsignedWord umiB;
    // int64  cbI;
    std::int32_t cbMatch;  //-1: no match, 0: exact, 1: 1 match with 1MM, >1: # of matches with 1MM
    std::int32_t umiCheck; // umi check status
    std::string cbMatchString;            // CB matches and qualities
    std::vector<UnsignedWord> cbMatchInd; // matches
    std::vector<std::uint32_t> cbReadCountExact;
    // map <uint32,uint32> cbReadCountMap;//count read per CB for no WL

    std::array<UnsignedWord, 256> qualHist;
    BarcodeStatistics stats;

    ReadBarcode(Context& context_in);
    void getCBandUMI(char** readSeq, char** readQual, UnsignedWord* readLen,
                     const std::string& readNameExtraIn, const std::uint32_t& readFilesIndex,
                     const char* readName);
    void addCounts(const ReadBarcode& rfIn);
    void addStats(const ReadBarcode& rfIn);
    void statsOut(std::ofstream& streamOut);
    void matchCBtoWL(std::string& cbSeq1, std::string& cbQual1, std::vector<UnsignedWord>& cbWL,
                     std::int32_t& cbMatch1, std::vector<UnsignedWord>& cbMatchInd1,
                     std::string& cbMatchString1);
    bool convertCheckUMI();
    void addStats(const std::int32_t cbMatch1);

protected:
    void readNameExtra();
};

} // namespace nustar::star2711b::solo::detail
