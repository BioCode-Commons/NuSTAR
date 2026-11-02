#pragma once

#include "../../stitching/paired_interface.hpp"

#include <memory>
#include "../../common/workspace.hpp"
#include "../../stitching/candidate.hpp"
#include "../../seed_search/interface.hpp"
#include "../../window/interface.hpp"
#include "../../stitching/config.hpp"
#include "../../stitching/score.hpp"
#include "../../common/compatibility_limits.hpp"

namespace nustar::star2711b {

class IPairedCandidateMapper : public IStageKernel {
public:
    virtual ~IPairedCandidateMapper() = default;
    virtual Result<std::unique_ptr<StageWorkspace>> createWorkspace(const WorkspaceBudget&) const = 0;
    // This is a joint paired search, not two finalized single-end alignments followed by pairing.
    virtual Status mapPair(const PairedMappingRead&, StageWorkspace&, CandidateSet& out) const = 0;
};

struct PairedCandidateMappingConfig {
    WindowConfig window;
    PairedStitchConfig stitch;
    ScoreConfig score;
    CompatibilityLimits limits;
};

Result<std::shared_ptr<const IPairedCandidateMapper>> composePairedCandidateMapper(
    ReferenceLease, std::shared_ptr<const ISeedFinder>,
    std::shared_ptr<const IWindowBuilder>, std::shared_ptr<const IPairedStitcher>,
    const PairedCandidateMappingConfig&);

} // namespace nustar::star2711b
