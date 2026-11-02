#pragma once

#include <memory>
#include "model.hpp"
#include "../common/workspace.hpp"
#include "../reads/model.hpp"
#include "../reference/view.hpp"
#include "../stitching/candidate.hpp"
#include "../selection/model.hpp"
#include "config.hpp"
#include "../selection/config.hpp"
#include "../stitching/config.hpp"
#include "../common/compatibility_limits.hpp"

namespace nustar::star2711b {

class ISingleSamProjector : public IStageKernel {
public:
    virtual ~ISingleSamProjector() = default;
    virtual Result<std::unique_ptr<StageWorkspace>> createWorkspace(
        const SamConfig&, const CompatibilityLimits&, const WorkspaceBudget&) const = 0;
    virtual Result<ProjectedOutput> projectSingle(const PreparedSingleRead&, const CandidateSet&,
        const SelectionResult&, const FilterDecision&, const ReferenceView&,
        const SamConfig&, const SelectionConfig&, const SingleStitchConfig&,
        StageWorkspace&) const = 0;
};

} // namespace nustar::star2711b
