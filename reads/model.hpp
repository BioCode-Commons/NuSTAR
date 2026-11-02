#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>
#include "raw.hpp"
#include "prepare_config.hpp"
#include <utility>

namespace nustar::star2711b {

struct ClipState {
    std::uint64_t fixed_5p, fixed_3p, adapter_5p, adapter_3p;
    std::uint64_t after_adapter_5p, after_adapter_3p;
    std::uint64_t total_5p, total_3p; // Original read coordinates, before strand conversion.
};
struct PreparedMate {
    std::uint64_t original_length, clipped_length;
    ClipState clipping;
    ReadRange combined_forward_interval;
    bool reverse_complemented_in_combined;
};
struct ReadBudget {
    std::uint64_t mismatch_limit; // min(abs limit, truncated fraction * sum(clipped mate lengths)).
    std::uint64_t legacy_Lread; // SE: L1. PE: L1 + 1 + L2, including one spacer.
    std::uint64_t legacy_original_pair_length;
    std::uint64_t biological_bases; // Sum of clipped mate lengths, excluding spacer.
};
// Normalized input for shared seed/window/selection/SJ kernels, without output records.
struct MappingRead {
    FragmentId fragment;
    std::vector<PreparedMate> mates; // Exactly one or two, supplied by the SE/PE preparer.
    std::array<std::vector<BaseCode>, 3> encoded; // Forward, complement, reverse complement.
    ReadBudget budget;
    bool revert_strand;
};
template<Layout L> class MappingInput {
public:
    // Wrap prepared geometry without revalidation; the preparer/PE merge owner sets it.
    // a merged single input inherits the original pair's mismatch budget.
    // Do not recompute inherited mismatch budgets for a PE-merged single mapping input.
    static MappingInput fromPrepared(MappingRead&&);
    const MappingRead& view() const noexcept;
private:
    explicit MappingInput(MappingRead&&);
    MappingRead read_;
};

template<Layout L>
MappingInput<L> MappingInput<L>::fromPrepared(MappingRead&& read) {
    return MappingInput<L>(std::move(read));
}

template<Layout L>
const MappingRead& MappingInput<L>::view() const noexcept {
    return read_;
}

template<Layout L>
MappingInput<L>::MappingInput(MappingRead&& read) : read_(std::move(read)) {}

using SingleMappingRead = MappingInput<Layout::SE>;
using PairedMappingRead = MappingInput<Layout::PE>;
struct PreparedSingleRead {
    RawSingleRead original;
    std::string output_name;
    std::string quality_for_output; // Converted FASTQ quality or STAR's synthetic FASTA quality.
    SingleMappingRead mapping;
};
struct PreparedReadPair {
    RawReadPair original;
    std::string output_name;
    std::array<std::string, 2> quality_for_output;
    PairedMappingRead mapping;
};

} // namespace nustar::star2711b
