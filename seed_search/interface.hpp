#pragma once

#include <memory>
#include <string_view>
#include "../common/workspace.hpp"
#include "delivery.hpp"
#include "plan.hpp"
#include "../index/backend.hpp"
#include "../reads/model.hpp"
#include "../common/compatibility_limits.hpp"

namespace nustar::star2711b {

class ISeedFinder : public IStageKernel {
public:
    virtual ~ISeedFinder() = default;
    virtual const SeedDeliveryContract& outputContract() const noexcept = 0;
    virtual Result<std::unique_ptr<StageWorkspace>> createWorkspace(
        const CompatibilityLimits&, const WorkspaceBudget&) const = 0;
    // Reference and SeedSearchPlan were frozen at bind(); they cannot change per call.
    // Caller supplies this finder's workspace and a prepared MappingRead.
    virtual Status find(const MappingRead&, StageWorkspace&, SeedBatch& out) const = 0;
};
class ISeedAlgorithmFactory {
public:
    virtual ~ISeedAlgorithmFactory() = default;
    virtual std::string_view algorithmId() const noexcept = 0;
    // Native index capabilities plus adapter behavior determine the delivered contract.
    virtual Result<std::shared_ptr<const ISeedFinder>> bind(
        SearchIndexLease, ReferenceLease, std::shared_ptr<const SeedSearchPlan>) const = 0;
};

} // namespace nustar::star2711b
