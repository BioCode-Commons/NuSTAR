#pragma once

#include "../common/workspace.hpp"
#include "candidate.hpp"
#include "../reads/model.hpp"
#include "../window/model.hpp"
#include "../reference/view.hpp"
#include "config.hpp"
#include "score.hpp"
#include "../common/compatibility_limits.hpp"

namespace nustar::star2711b {

class ISingleStitcher : public IStageKernel {
public:
    virtual ~ISingleStitcher() = default;
    virtual Result<std::unique_ptr<StageWorkspace>> createWorkspace(
        const SingleStitchConfig&, const ScoreConfig&, const ReferenceView&,
        const CompatibilityLimits&, const WorkspaceBudget&) const = 0;
    // Inputs/reference/config/workspace must belong to the same bound composition.
    virtual Status stitchSingle(const SingleMappingRead&, const WindowBatch&, const ReferenceView&,
        const SingleStitchConfig&, const ScoreConfig&, StageWorkspace&, CandidateSet& out) const = 0;
};

} // namespace nustar::star2711b
