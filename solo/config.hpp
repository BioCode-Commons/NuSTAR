#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace nustar::star2711b::solo {

struct Config {
    std::string type = "None";
    unsigned barcode_mate = 0;
    bool bam_input = false, coordinate_bam = false;
    std::vector<std::string> output_attributes;
    std::array<std::array<std::uint64_t, 2>, 2> clipping{};
    std::vector<std::string> read_groups;
    std::vector<std::string> read_group_lines;
    // STAR option names and tokens, in input order.
    std::vector<std::string> arguments;
};

} // namespace nustar::star2711b::solo
