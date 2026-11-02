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

class IPairedStitcher : public IStageKernel {
public:
    virtual ~IPairedStitcher() = default;
    virtual Result<std::unique_ptr<StageWorkspace>> createWorkspace(
        const PairedStitchConfig&, const ScoreConfig&, const ReferenceView&,
        const CompatibilityLimits&, const WorkspaceBudget&) const = 0;
    // Joint search must preserve mate gaps, overlap consistency, and per-mate score retention.
    // Inputs/reference/config/workspace must belong to the same bound composition.
    virtual Status stitchPair(const PairedMappingRead&, const WindowBatch&, const ReferenceView&,
        const PairedStitchConfig&, const ScoreConfig&, StageWorkspace&, CandidateSet& out) const = 0;
};

} // namespace nustar::star2711b
