#pragma once

#include <optional>
#include "delivery.hpp"
#include "../index/backend.hpp"
#include "config.hpp"

namespace nustar::star2711b {

struct SeedSearchPlan {
    SeedDeliveryContract delivery;
    // A versioned semantic search universe, separate from physical index coverage and tuning.
    std::optional<BackendOptions> equivalent_search_domain;
    SeedConfig compatibility_policy;
    BackendOptions algorithm_options;
};

} // namespace nustar::star2711b
