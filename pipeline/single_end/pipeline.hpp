#pragma once

#include <memory>
#include "mapping_interface.hpp"
#include "../../selection/interface.hpp"
#include "../../sam/single_projector.hpp"
#include "result.hpp"
#include "../effects.hpp"
#include "../../splice_junction/config.hpp"

namespace nustar::star2711b {

struct SingleKernelBundle {
    std::shared_ptr<const ISingleCandidateMapper> mapper;
    std::shared_ptr<const IAlignmentSelector> selection;
    std::shared_ptr<const ISingleSamProjector> sam;
};

struct SinglePipelineConfig {
    SelectionConfig select;
    ReadFilterConfig filter;
    SjOutputConfig sj;
    SamConfig sam;
    SingleStitchConfig stitch;
    CompatibilityLimits limits;
};

class ISingleReadPipeline : public IStageKernel {
public:
    virtual ~ISingleReadPipeline() = default;
    virtual Result<std::unique_ptr<StageWorkspace>> createWorkspace(const WorkspaceBudget&) const = 0;
    virtual Result<FragmentEffects> computeSingle(
        PreparedSingleRead&&, StageWorkspace&, LegacyRandomState&) const = 0;
    // Diagnostic path captures owned read/candidate/selection values from the same stages.
    // Validation stage wrappers separately capture seed and window snapshots.
    virtual Result<SingleReadComputation> traceSingle(
        PreparedSingleRead&&, StageWorkspace&, LegacyRandomState&) const = 0;
};

// Independent orchestration implementations; neither factory takes the other pipeline.
Result<std::shared_ptr<const ISingleReadPipeline>> composeSinglePipeline(
    ReferenceLease, SingleKernelBundle, const SinglePipelineConfig&);

} // namespace nustar::star2711b
