#pragma once

#include <memory>
#include "../reads/input_config.hpp"
#include "../reads/prepare_config.hpp"
#include "backend_names.hpp"
#include "resources.hpp"
#include "../seed_search/plan.hpp"
#include "../pipeline/paired_end/pipeline.hpp"

namespace nustar::star2711b {

struct PairedRunPlan {
    std::shared_ptr<const RunResources> resources;
    InputSpec input; // Validated PE transport configuration.
    PairedPrepareConfig prepare;
    PairedBackendNames resolved_backends;
    std::shared_ptr<const SeedSearchPlan> merged_seed_search;
    std::shared_ptr<const IPairedReadPipeline> pipeline;
};

} // namespace nustar::star2711b
