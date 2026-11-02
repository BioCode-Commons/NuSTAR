#pragma once

#include <memory>
#include "../common/workspace.hpp"
#include "model.hpp"
#include "../seed_search/delivery.hpp"
#include "../reads/model.hpp"
#include "../reference/view.hpp"
#include "config.hpp"
#include "../common/compatibility_limits.hpp"

namespace nustar::star2711b {

class IWindowBuilder : public IStageKernel {
public:
    virtual ~IWindowBuilder() = default;
    virtual Status validateSeedContract(const SeedDeliveryContract&, const WindowConfig&) const = 0;
    virtual Result<std::unique_ptr<StageWorkspace>> createWorkspace(
        const WindowConfig&, const ReferenceView&, const CompatibilityLimits&,
        const WorkspaceBudget&) const = 0;
    // Composition validates the seed contract once. Inputs, reference, config,
    // and workspace must belong to that composition and the same fragment.
    virtual Status build(const MappingRead&, const SeedBatch&, const ReferenceView&,
        const WindowConfig&, StageWorkspace&, WindowBatch& out) const = 0;
};

} // namespace nustar::star2711b
