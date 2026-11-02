#include "planner.hpp"
#include "../common/unavailable.hpp"

#include <utility>

namespace nustar::star2711b {

SinglePassPlanner::SinglePassPlanner(std::shared_ptr<const IExecutionAssembler> assembler)
    : assembler_(std::move(assembler)) {}

Result<PlanningDecision> SinglePassPlanner::prepare(
    const RunRequest&, IReferenceService&, IIndexBackend&) {
    return unavailable(Stage::Plan, assembler_ ? "run planner" : "execution assembler");
}

} // namespace nustar::star2711b
