#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include "../single_end/mapping_interface.hpp"
#include "overlap_state.hpp"
#include "overlap_config.hpp"

namespace nustar::star2711b {

// Compositions freeze dependencies and never resolve backends during mapping.
struct MergedPairMapping {
    SingleMappingRead merged;
    std::array<ReadOffset, 2> mate_starts_in_merged;
    std::uint64_t overlapping_bases;
};
class IPairedOverlapRefiner : public IStageKernel {
public:
    virtual ~IPairedOverlapRefiner() = default;
    virtual Result<std::unique_ptr<StageWorkspace>> createWorkspace(
        const OverlapConfig&, const ISingleCandidateMapper&, const CompatibilityLimits&,
        const WorkspaceBudget&) const = 0;
    // PE owns merge/projection/rescoring. The helper maps only a SingleMappingRead.
    // Its workspace is separate from main mapping; it cannot invoke a single-end pipeline.
    virtual Status refinePair(const PairedMappingRead&, const ReferenceView&,
        const OverlapConfig&, const ScoreConfig&, const ISingleCandidateMapper&,
        StageWorkspace&, CandidateSet& inout, OverlapState& out) const = 0;
};

Result<std::shared_ptr<const IPairedOverlapRefiner>> makeStarPairedOverlapRefiner(ReferenceLease);

} // namespace nustar::star2711b
