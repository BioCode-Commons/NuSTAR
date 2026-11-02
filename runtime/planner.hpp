#pragma once

// Reserved framework contract; see README for current implementation status.

#include <memory>
#include "plan.hpp"

namespace nustar::star2711b {

class SinglePassPlanner {
public:
    explicit SinglePassPlanner(std::shared_ptr<const IExecutionAssembler>);
    // May choose a whole-run legacy path after index-derived options are resolved.
    // No formal alignment output is opened before this decision.
    Result<PlanningDecision> prepare(const RunRequest&, IReferenceService&,
        IIndexBackend&);
private:
    std::shared_ptr<const IExecutionAssembler> assembler_; // Explicit registry/composition owner.
};

} // namespace nustar::star2711b
