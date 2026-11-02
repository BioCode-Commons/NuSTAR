#pragma once

#include <cstdint>
#include <filesystem>
#include "../common/types.hpp"

namespace nustar::star2711b {

struct RuntimeConfig {
    std::uint32_t runThreadN, runRNGseed;
    std::uint64_t input_buffer_bytes, output_buffer_bytes;
    std::uint64_t queued_output_bytes_max;
    std::filesystem::path output_prefix, temporary_directory;
};

} // namespace nustar::star2711b
