#pragma once

#include "model.hpp"

namespace nustar::star2711b {

// Takes ownership of transport strings; default product FASTQ preparation; configurable clipping remains a separate contract.
PreparedSingleRead prepareSingleFastq(std::string header, std::string sequence,
    std::string quality, std::uint64_t ordinal,
    std::array<std::uint64_t,2> clipping={}, std::uint32_t file_group=0);
PreparedReadPair preparePairedFastq(std::string header1, std::string sequence1,
    std::string quality1, std::string header2, std::string sequence2,
    std::string quality2, std::uint64_t ordinal,
    std::array<std::array<std::uint64_t,2>,2> clipping={}, std::uint32_t file_group=0);

} // namespace nustar::star2711b
