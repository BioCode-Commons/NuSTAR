#pragma once

#include <memory>
#include "mapping_interface.hpp"
#include "../single_end/mapping_interface.hpp"
#include "overlap.hpp"
#include "../../selection/interface.hpp"
#include "../../sam/paired_projector.hpp"
#include "result.hpp"
#include "../effects.hpp"
#include "../../splice_junction/config.hpp"

namespace nustar::star2711b {

struct PairedKernelBundle {
    std::shared_ptr<const IPairedCandidateMapper> mapper;
    std::shared_ptr<const ISingleCandidateMapper> merged_single_mapper;
    std::shared_ptr<const IPairedOverlapRefiner> overlap;
    std::shared_ptr<const IAlignmentSelector> selection;
    std::shared_ptr<const IPairedSamProjector> sam;
};

struct PairedPipelineConfig {
    SelectionConfig select;
    ReadFilterConfig filter;
    SjOutputConfig sj;
    PairedSamConfig sam;
    PairedStitchConfig stitch;
    OverlapConfig overlap;
    ScoreConfig score;
    CompatibilityLimits limits;
};

class IPairedReadPipeline : public IStageKernel {
public:
    virtual ~IPairedReadPipeline() = default;
    virtual Result<std::unique_ptr<StageWorkspace>> createWorkspace(const WorkspaceBudget&) const = 0;
    virtual Result<FragmentEffects> computePair(
        PreparedReadPair&&, StageWorkspace&, LegacyRandomState&) const = 0;
    virtual Result<PairedReadComputation> tracePair(
        PreparedReadPair&&, StageWorkspace&, LegacyRandomState&) const = 0;
};

Result<std::shared_ptr<const IPairedReadPipeline>> composePairedPipeline(
    ReferenceLease, PairedKernelBundle, const PairedPipelineConfig&);

} // namespace nustar::star2711b
