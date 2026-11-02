#pragma once

#include "../counting/record_types.hpp"
#include "../support/types.hpp"

namespace nustar::star2711b::solo::detail {

class Options;

class BarcodeLayout { // complex barcodes
    //    ParametersSolo *pSolo;
public:
    std::uint32_t mate; // 0 or 1
    // array: 0=start, 1=end
    std::int32_t anchorType[2]; // 0=read start; 1=read end; 2=adapter start; 3=adapter end
    std::int32_t anchorDist[2]; // distance to anchor

    std::uint32_t adapterLength; // length of the adapter

    std::vector<std::vector<BarcodeCode>> wl; // whitelists, one for each length
    std::vector<std::vector<BarcodeCode>>
        wlEd; // edited whitelists (i.e. including mismatches and indels)
    std::vector<std::vector<std::uint32_t>> wlEdInd; // index for wlEd in the unedited wl

    UnsignedWord wlFactor; // factor and modulo for converting each whitelist index into global
                           // index
    std::vector<std::uint32_t> wlAdd; // additive for each length
    std::uint32_t minLen;             // minimum length for this barcode
    std::uint32_t totalSize;          // total size of all whitelists

    UnsignedWord iLen, iCB; // indexes, used in ParametersSolo.cpp

    bool extractBarcode(std::string& seqIn, std::string& qualIn, const std::uint32_t aStart,
                        std::string& bSeq, std::string& bQual);
    void sortWhiteList(Options* options);
    void extractPositionsFromString(std::string& strIn);

    // SoloBarcode(ParametersSolo *pSolo) : pSolo(pSolo) {};
};

} // namespace nustar::star2711b::solo::detail
