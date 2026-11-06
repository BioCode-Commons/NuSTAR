#pragma once

#include "../index/backend.hpp"
#include "../bam/config.hpp"
#include "../reads/input_config.hpp"
#include "../reference/genome.hpp"
#include "../solo/config.hpp"
#include "../pipeline/paired_end/overlap_config.hpp"
#include <limits>

namespace nustar::star2711b {

struct AlignmentResources {
    std::shared_ptr<const GenomeReferenceData> reference;
    SearchIndexLease index;
};

struct RunControl {
    OverlapConfig overlap{0,0.01};
    AlignmentResources resources;
    std::vector<InputGroup> input_groups;
    std::array<std::array<std::uint64_t, 2>, 2> clipping{};
    std::shared_ptr<const solo::Config> solo_config;
    std::optional<std::filesystem::path> barcode_fastq;
    bool sam_enabled = true;
    bool bam_input = false;
    std::vector<std::string> input_attributes{"All"};
    BamConfig bam;
    std::vector<std::string> output_attributes;
    bool twopass = false;
    bool save_inserted_index = false;
    std::uint64_t pass1_reads = std::numeric_limits<std::uint64_t>::max();
    std::uint64_t read_limit = std::numeric_limits<std::uint64_t>::max();
    std::uint64_t sjdb_overhang = 100;
    std::uint64_t junction_limit = 1000000;
};

} // namespace nustar::star2711b
