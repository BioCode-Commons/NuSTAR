#pragma once

#include <memory>
#include "../reference/view.hpp"
#include "../index/backend.hpp"
#include "../seed_search/plan.hpp"

namespace nustar::star2711b {

struct ResolvedRunConfig;

struct RunResources {
    ReferenceLease reference;
    SearchIndexLease search_index;
    std::shared_ptr<const SeedSearchPlan> seed_search; // Shared with the bound finder.
    std::shared_ptr<const ResolvedRunConfig> config;
};

} // namespace nustar::star2711b
