#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>
#include "../common/types.hpp"

namespace nustar::star2711b {

inline bool isSoloAttribute(std::string_view tag) {
    return tag=="CR" || tag=="CY" || tag=="UR" || tag=="UY" || tag=="CB" || tag=="UB" ||
        tag=="GX" || tag=="GN" || tag=="sM" || tag=="sS" || tag=="sQ" || tag=="sF";
}
struct SamConfig {
    bool enabled;
    std::string outSAMmode, outSAMorder, outSAMreadID;
    std::vector<std::string> outSAMattributes; // Expanded tags in exact output order.
    std::uint64_t outSAMattrIHstart;
    std::uint8_t outSAMmapqUnique;
    std::uint16_t outSAMflagOR, outSAMflagAND;
    bool unmapped_within, unmapped_fastx;
    std::vector<std::string> read_group_lines;
    std::vector<std::string> read_group_by_file;
    std::string header_hd, header_pg; // User-specified HD/PG lines, not the generated STAR PG line.
    std::string command_line_full, command_line_user; // Generated @PG CL and @CO respectively.
    std::string header_extra; // Exact resolved extra header bytes; empty for ordinary Fastx input.
    std::vector<std::string> extra_reference_lines; // Exact extraReferences.txt @SQ lines, in file order.
    std::vector<std::string> header_comments; // User comment-file lines; excludes generated @CO.
    std::uint64_t limitOutSAMoneReadBytes;
};
struct PairedSamConfig {
    SamConfig common;
    bool keep_pairs;
};

} // namespace nustar::star2711b
