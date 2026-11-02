#pragma once

#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <variant>
#include <vector>
#include "run_config.hpp"
#include "../reads/input_config.hpp"
#include "reference_service.hpp"
#include "../seed_search/plan.hpp"
#include "../index/backend.hpp"
#include "backend_names.hpp"
#include "report.hpp"
#include "single_plan.hpp"
#include "paired_plan.hpp"

namespace nustar::star2711b {

struct RunRequest {
    InputSpec input;
    ReferenceRequest reference;
    std::optional<ReferenceAnnotationRequest> annotation_update;
    SearchIndexRequest search_index;
    SeedSearchPlan seed_search;
    // Null inherits the main seed plan during planning; the resolved auxiliary binding is recorded.
    std::optional<SeedSearchPlan> paired_merged_seed_search;
    std::vector<std::string> original_arguments;
    std::vector<std::filesystem::path> parameter_files;
    BackendNames requested_backends;
};
using RunPlan = std::variant<SingleRunPlan, PairedRunPlan>;
struct LegacyRunPlan {
    RunRequest request; // The complete original invocation is passed to the legacy entry.
    std::vector<OptionDecision> option_decisions; // Includes the reason for whole-run fallback.
};
using PlanningDecision = std::variant<RunPlan, LegacyRunPlan>;
// Runtime assembly may inspect full config; algorithm kernels and workspaces may not.
class IExecutionAssembler {
public:
    virtual ~IExecutionAssembler() = default;
    virtual Result<SingleRunPlan> assembleSingle(std::shared_ptr<const RunResources>,
        const InputSpec&, const SingleBackendNames&) const = 0;
    virtual Result<PairedRunPlan> assemblePaired(std::shared_ptr<const RunResources>,
        const InputSpec&, const PairedBackendNames&,
        std::shared_ptr<const SeedSearchPlan> merged_seed_search) const = 0;
};

} // namespace nustar::star2711b
