#pragma once
#include "../support/types.hpp"

#include <array>
#include <unordered_map>

namespace nustar::star2711b::solo::detail {

typedef struct {
    UnsignedWord cb;
    std::uint32_t umi;
} ReadCorrection;

typedef struct {
    std::uint32_t tr;
    std::uint8_t type;
} TranscriptType;

typedef std::uint32_t UmiCode;
typedef UnsignedWord BarcodeCode;
typedef std::uint32_t ReadCount;

inline constexpr std::uint32_t uintUMIbits = 32;
inline constexpr std::uint32_t velocytoTypeGeneBits = 4;
inline constexpr std::uint32_t velocytoTypeGeneBitShift = 28;
inline constexpr std::uint32_t geneMultMark = (((std::uint32_t)1) << 31);

class ReadFlags {
public:
    bool yes = false;
    typedef std::uint32_t typeFlag;
    typeFlag flag = 0;
    enum : std::uint32_t {
        cbMatch,
        cbPerfect,
        cbMMunique,
        cbMMmultiple,
        genomeU,
        genomeM,
        featureU,
        featureM,
        exonic,
        intronic,
        exonicAS,
        intronicAS,
        mito,
        countedU,
        countedM,
        nBits
    };
    const std::vector<std::string> statNames = {
        "cbMatch",  "cbPerfect",  "cbMMunique", "cbMMmultiple", "genomeU",
        "genomeM",  "featureU",   "featureM",   "exonic",       "intronic",
        "exonicAS", "intronicAS", "mito",       "countedU",     "countedM"};

    

    std::unordered_map<BarcodeCode, std::array<UnsignedWord, nBits>> flagCounts;
    std::array<UnsignedWord, nBits> flagCountsNoCB = {};

    void setBit(std::uint32_t ibit) {
        flag |= ((typeFlag)1) << ibit;
    };

    typeFlag checkBit(std::uint32_t ibit) {
        return (flag >> ibit) & ((typeFlag)1);
    };

    void countsAdd(BarcodeCode cb) { // adds flag bits to the count for a given cb
        auto cbInserted = flagCounts.insert({cb, {}});
        for (std::uint32_t ibit = 0; ibit < nBits; ibit++) {
            (*cbInserted.first).second[ibit] += (UnsignedWord)checkBit(ibit);
        };
    };

    void countsAddNoCB() { // adds flag bits to the count for no-CB reads
        for (std::uint32_t ibit = 0; ibit < nBits; ibit++)
            flagCountsNoCB[ibit] += (UnsignedWord)checkBit(ibit);
    };

    void countsAddNoCBarray(std::array<UnsignedWord, nBits>& arrIn) {
        for (std::uint32_t ibit = 0; ibit < nBits; ibit++)
            flagCountsNoCB[ibit] += arrIn[ibit];
    };
};

} // namespace nustar::star2711b::solo::detail
