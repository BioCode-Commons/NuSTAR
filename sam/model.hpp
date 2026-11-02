#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>
#include "../common/types.hpp"

namespace nustar::star2711b {

struct CigarOperation { std::uint64_t length; char operation; };
struct SamAuxField {
    std::array<char, 2> tag;
    char type; // SAM type letter; numeric array subtype stays inside encoded_value.
    std::string encoded_value; // Exact STAR text; no lossy parsing or canonicalization.
};
struct SamRecord {
    FragmentId fragment;
    MateId mate;
    std::uint64_t alignment_ordinal;
    std::string qname;
    std::uint16_t flag;
    std::string rname;
    std::uint64_t position_1based; // Zero means unavailable, as required by SAM.
    std::uint8_t mapq;
    std::vector<CigarOperation> cigar; // Empty encodes '*'.
    std::string rnext; // Preserve '*', '=', or the exact reference name.
    std::uint64_t mate_position_1based;
    std::int64_t template_length;
    std::string sequence, quality; // Already oriented; '*' is preserved literally.
    std::vector<SamAuxField> tags;
};
struct UnmappedFastxRecord { MateId mate; std::string encoded_record; };
struct OutputDisposition {
    std::array<bool, 2> mate_mapped;
    int legacy_unmap_type_after_output; // Includes one-mate-mapped PE state 4.
    bool count_as_unmapped_all;
};
struct ProjectedOutput {
    std::vector<SamRecord> records;
    std::vector<UnmappedFastxRecord> unmapped;
    OutputDisposition disposition;
};

} // namespace nustar::star2711b
