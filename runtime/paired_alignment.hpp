#pragma once

#include "run_control.hpp"
#include <cstdint>
#include <filesystem>
#include <string>
#include <optional>
#include "../common/types.hpp"

namespace nustar::star2711b {

struct PairedSamOptions {
    std::filesystem::path genome_dir;
    std::filesystem::path mate1_fastq, mate2_fastq;
    std::filesystem::path output_sam;
    std::string command_line;
    std::uint32_t threads = 1;
    std::string io_backend = "mmap";
    // Empty preserves the legacy explicit-SAM-path API's auxiliary filenames.
    std::optional<std::string> outFileNamePrefix;
    bool outSAMunmappedWithin = true;
    std::string outSAMorder = "PairedKeepInputOrder";
    RunControl run;
};

struct PairedSamReport {
    std::uint64_t pairs = 0;
    std::uint64_t records = 0;
    double initialization_seconds = 0, mapping_wall_seconds = 0;
};

Result<PairedSamReport> runPairedSam(const PairedSamOptions&);

} // namespace nustar::star2711b
