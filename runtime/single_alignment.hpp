#pragma once

#include "run_control.hpp"
#include <cstdint>
#include <filesystem>
#include <string>
#include <optional>
#include "../common/types.hpp"

namespace nustar::star2711b {

// Narrow, validated single-end FASTQ-to-SAM/SJ/statistics path. The full
// runner with broader input and STAR option support is a separate contract.
struct SingleSamOptions {
    std::filesystem::path genome_dir;
    std::filesystem::path reads_fastq;
    std::filesystem::path output_sam;
    std::string command_line;
    bool profile_stages = false;
    std::uint32_t threads = 1;
    std::string io_backend = "mmap";
    // Empty preserves the legacy explicit-SAM-path API's auxiliary filenames.
    std::optional<std::string> outFileNamePrefix;
    bool outSAMunmappedWithin = true;
    std::string outSAMorder = "PairedKeepInputOrder";
    RunControl run;
};

struct SingleSamReport {
    std::uint64_t reads = 0;
    std::uint64_t records = 0;
    double seed_seconds = 0;
    double window_seconds = 0;
    double stitch_seconds = 0;
    double selection_seconds = 0;
    double sam_seconds = 0;
    // Initialization/mapping are wall times; waits are sums over workers.
    double reference_load_seconds = 0;
    double index_load_seconds = 0;
    double seed_bind_seconds = 0;
    double pipeline_setup_seconds = 0;
    double mapping_wall_seconds = 0;
    double input_wait_seconds = 0;
    double commit_wait_seconds = 0;
    double prepare_seconds = 0;
    double input_read_seconds = 0;
};

Result<SingleSamReport> runSingleSam(const SingleSamOptions&);

} // namespace nustar::star2711b
