#pragma once

#include "../../stitching/single_interface.hpp"

#include <memory>
#include "../../common/workspace.hpp"
#include "../../stitching/candidate.hpp"
#include "../../seed_search/interface.hpp"
#include "../../window/interface.hpp"
#include "../../stitching/config.hpp"
#include "../../stitching/score.hpp"
#include "../../common/compatibility_limits.hpp"

namespace nustar::star2711b {

class ISingleCandidateMapper : public IStageKernel {
public:
    virtual ~ISingleCandidateMapper() = default;
    // Mapping only; excludes selection, overlap, output, statistics, and RNG consumption.
    virtual Result<std::unique_ptr<StageWorkspace>> createWorkspace(const WorkspaceBudget&) const = 0;
    virtual Status mapSingle(const SingleMappingRead&, StageWorkspace&, CandidateSet& out) const = 0;
};

struct SingleCandidateMappingConfig {
    WindowConfig window;
    SingleStitchConfig stitch;
    ScoreConfig score;
    CompatibilityLimits limits;
};

Result<std::shared_ptr<const ISingleCandidateMapper>> composeSingleCandidateMapper(
    ReferenceLease, std::shared_ptr<const ISeedFinder>,
    std::shared_ptr<const IWindowBuilder>, std::shared_ptr<const ISingleStitcher>,
    const SingleCandidateMappingConfig&);

} // namespace nustar::star2711b
