#pragma once

#include <memory>
#include "../reads/input_config.hpp"
#include "../reads/prepare_config.hpp"
#include "backend_names.hpp"
#include "resources.hpp"
#include "../pipeline/single_end/pipeline.hpp"

namespace nustar::star2711b {

struct SingleRunPlan {
    std::shared_ptr<const RunResources> resources;
    InputSpec input; // Validated SE transport configuration.
    SinglePrepareConfig prepare;
    SingleBackendNames resolved_backends;
    std::shared_ptr<const ISingleReadPipeline> pipeline;
};

} // namespace nustar::star2711b
