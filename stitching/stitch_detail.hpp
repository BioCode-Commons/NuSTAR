#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace nustar::star2711b::stitch_detail {

using uint = unsigned long long;
using intScore = int;
using namespace std;

constexpr unsigned MAX_N_EXONS = 20;
constexpr unsigned MAX_SJ_REPEAT_SEARCH = 255;
constexpr uint DEF_readSeqLengthMax = 650;
constexpr unsigned char MARK_FRAG_SPACER_BASE = 11;
constexpr int scoreMatch = 1;
constexpr int EXIT_CODE_INPUT_FILES = 1;

#define EXTEND_ORDER 1
constexpr std::size_t EX_R = 0, EX_G = 1, EX_L = 2, EX_iFrag = 3, EX_sjA = 4;
constexpr std::size_t WA_Length = 0, WA_rStart = 1, WA_gStart = 2;
constexpr std::size_t WA_Nrep = 3, WA_Anchor = 4, WA_iFrag = 5, WA_sjA = 6;

using uiWA = std::array<uint, 7>;

struct Genome {
    const char* G = nullptr;
    uint nGenome = 0;
    std::vector<uint> chrStart, chrLength;
    std::vector<uint> sjdbStart, sjdbEnd;
    std::vector<unsigned char> sjdbMotif, sjdbShiftLeft, sjdbShiftRight, sjdbStrand;
    uint sjdbN = 0;
};

struct Parameters {
    struct { bool yes = true; } alignSoftClipAtReferenceEnds;
    struct { bool flushRight = false; } alignInsertionFlush;
    struct { std::array<std::array<bool, 2>, 2> ext{}; std::string in; } alignEndsType;
    struct { uint nBasesMax = 0; } alignEndsProtrude;
    struct { int type = 0; } outSAMstrandField;
    struct { int segmentMin = 0; } pCh;
    struct { int sjdbScore = 0; } pGe;

    uint alignIntronMin = 21, alignIntronMax = 0;
    uint alignSJoverhangMin = 5, alignSJDBoverhangMin = 3;
    std::array<int, 4> alignSJstitchMismatchNmax{};
    uint alignSplicedMateMapLmin = 0;
    double alignSplicedMateMapLminOverLmate = 0;
    uint alignTranscriptsPerWindowNmax = 10;
    uint alignMatesGapMax = 0;
    std::string outFilterIntronMotifs = "None", outFilterIntronStrands = "None";
    double outFilterMismatchNoverLmax = 0.3;
    int outFilterMultimapScoreRange = 1;
    int outFilterBySJoutStage = 0;
    uint* sjNovelStart = nullptr;
    uint* sjNovelEnd = nullptr;
    int sjNovelN = 0;
    int scoreGap = 0, scoreGapNoncan = 0, scoreGapGCAG = 0, scoreGapATAC = 0;
    int scoreDelOpen = 0, scoreDelBase = 0, scoreInsOpen = 0, scoreInsBase = 0;
    int scoreStitchSJshift = 0;
    double scoreGenomicLengthLog2scale = 0;
};

struct Transcript {
    uint exons[MAX_N_EXONS][5]{};
    uint shiftSJ[MAX_N_EXONS][2]{};
    int canonSJ[MAX_N_EXONS]{};
    unsigned char sjAnnot[MAX_N_EXONS]{}, sjStr[MAX_N_EXONS]{};
    uint intronMotifs[3]{};
    unsigned char sjMotifStrand = 0;
    bool sjYes = false;
    uint nExons = 0;
    int iFrag = 0;
    uint rStart = 0, roStart = 0, rLength = 0, gStart = 0, gLength = 0;
    uint Chr = 0, Str = 0, roStr = 0;
    uint nMatch = 0, nMM = 0, mappedLength = 0, extendL = 0;
    intScore maxScore = 0;
    uint nGap = 0, lGap = 0, nDel = 0, lDel = 0, nIns = 0, lIns = 0;
    uint nUnique = 0, nAnchor = 0;

    void reset() {
        extendL = rStart = roStart = rLength = gStart = gLength = 0;
        maxScore = 0;
        nMatch = nMM = nGap = lGap = nDel = lDel = nIns = lIns = 0;
        nUnique = nAnchor = 0;
    }
    void add(const Transcript* other) {
        maxScore += other->maxScore;
        nMatch += other->nMatch; nMM += other->nMM;
        nGap += other->nGap; lGap += other->lGap;
        nDel += other->nDel; lDel += other->lDel;
        nIns += other->nIns; lIns += other->lIns;
        nUnique += other->nUnique;
    }
    int variationAdjust(const Genome&, const char*) const { return 0; }
};

struct ReadAlign {
    uint outFilterMismatchNmaxTotal = 0;
    std::array<uint, 2> readLength{};
    std::array<int, 2> maxScoreMate{};
    uint iRead = 0;
};

inline int binarySearch2(uint start, uint end, const uint* starts,
                         const uint* ends, int count) {
    int low = 0, high = count;
    while (low < high) {
        const int middle = low + (high-low)/2;
        if (starts[middle] < start ||
            (starts[middle] == start && ends[middle] < end)) low = middle+1;
        else high = middle;
    }
    return low < count && starts[low] == start && ends[low] == end ? low : -1;
}
inline int binarySearch2(uint start, uint end, const std::vector<uint>& starts,
                         const std::vector<uint>& ends, int count) {
    return binarySearch2(start, end, starts.data(), ends.data(), count);
}
inline void exitWithError(const std::string& message) {
    throw std::invalid_argument(message);
}
uint blocksOverlap(const Transcript& first, const Transcript& second);
bool extendAlign(const char* read, const char* genome, uint read_start,
                 uint genome_start, int read_step, int genome_step, uint length,
                 uint previous_length, uint previous_mismatches,
                 uint maximum_mismatches, double mismatch_fraction,
                 bool extend_to_end, Transcript* result);
intScore stitchAlignToTranscript(uint rAend, uint gAend, uint rBstart,
    uint gBstart, uint length, uint fragment, uint sjdb_id, Parameters& parameters,
    const char* read, const Genome& genome, Transcript* transcript,
    uint mismatch_limit);
void stitchWindowAligns(uint index, uint count, int score, bool included[],
    uint read_end, uint genome_end, Transcript transcript, uint read_length,
    uiWA* alignments, const char* read, const Genome& genome, Parameters& parameters,
    Transcript** results, uint* result_count, ReadAlign* aligner);

} // namespace nustar::star2711b::stitch_detail
