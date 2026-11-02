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

class IPairedSamProjector : public IStageKernel {
public:
    virtual ~IPairedSamProjector() = default;
    virtual Result<std::unique_ptr<StageWorkspace>> createWorkspace(
        const PairedSamConfig&, const CompatibilityLimits&, const WorkspaceBudget&) const = 0;
    // Owns mate flags, RNEXT/PNEXT/TLEN, KeepPairs, and one-mate-mapped disposition.
    virtual Result<ProjectedOutput> projectPair(const PreparedReadPair&, const CandidateSet&,
        const SelectionResult&, const FilterDecision&, const ReferenceView&,
        const PairedSamConfig&, const SelectionConfig&, const PairedStitchConfig&,
        StageWorkspace&) const = 0;
};

} // namespace nustar::star2711b
