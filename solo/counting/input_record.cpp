#include "input_record.hpp"
#include "feature_recorder.hpp"
#include "../support/binary_search.hpp"

namespace nustar::star2711b::solo::detail {

bool soloInputFeatureUMI(std::fstream* strIn, std::int32_t featureType, bool readInfoYes,
                         std::array<std::vector<UnsignedWord>, 2>& sjAll, UnsignedWord& iread,
                         std::int32_t& cbmatch, std::uint32_t& feature, UnsignedWord& umi,
                         std::vector<std::uint32_t>& featVecU32, ReadFlags& readFlagCounts) {
    if (!(*strIn >> umi)) // end of file
        return false;

    if (readInfoYes) {
        *strIn >> iread >> readFlagCounts.flag;
    };

    switch (featureType) {
    case FeatureTypes::Gene:
    case FeatureTypes::GeneFull:
    case FeatureTypes::GeneFull_Ex50pAS:
    case FeatureTypes::GeneFull_ExonOverIntron:
        *strIn >> feature;
        break;

    case FeatureTypes::SJ:
        UnsignedWord sj[2];
        *strIn >> sj[0] >> sj[1];
        feature = (std::uint32_t)binarySearch2(sj[0], sj[1], sjAll[0].data(), sjAll[1].data(),
                                               sjAll[0].size());
        break;

    case FeatureTypes::Transcript3p:
        feature = 0;
        std::uint32_t ntr, in1;
        *strIn >> ntr;
        featVecU32.resize(2 * ntr);
        for (std::uint32_t ii = 0; ii < 2 * ntr; ii++) {
            *strIn >> in1;
            featVecU32[ii] = in1;
        };
        break;
    };

    *strIn >> cbmatch;

    return true;
};

} // namespace nustar::star2711b::solo::detail
