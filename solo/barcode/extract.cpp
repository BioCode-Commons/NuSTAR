#include "layout.hpp"
#include <array>

namespace nustar::star2711b::solo::detail {

bool BarcodeLayout::extractBarcode(
    std::string& seqIn, std::string& qualIn, const std::uint32_t adapterStart, std::string& bSeq,
    std::string& bQual) { // input: sequence seqIn, adapter start adapterStart
                          // output: start position of the barcode
    std::array<std::int32_t, 2> pos = {0, 0};
    for (std::uint32_t ii = 0; ii < 2; ii++) {
        switch (anchorType[ii]) { // this calculates the position of the anchor base
        case 0:                   // read start
            pos[ii] = 0;
            break;
        case 1: // read end
            pos[ii] = seqIn.size() - 1;
            break;
        case 2: // adapter start
            pos[ii] = (std::int32_t)adapterStart;
            break;
        case 3: // adapter end
            pos[ii] = (std::int32_t)adapterStart + adapterLength - 1;
            break;
        };
        pos[ii] += anchorDist[ii];
    };

    bSeq = "";
    bQual = "";
    if (pos[0] < 0 || pos[1] > (std::int32_t)seqIn.size() || pos[0] > pos[1]) // something went
                                                                              // wrong
        return false;

    bSeq = seqIn.substr(pos[0], pos[1] - pos[0] + 1);
    bQual = qualIn.substr(pos[0], pos[1] - pos[0] + 1);
    return true;
};

} // namespace nustar::star2711b::solo::detail
