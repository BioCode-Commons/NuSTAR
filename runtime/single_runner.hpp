#pragma once

// Reserved framework contract; see README for current implementation status.

#include <memory>
#include "single_plan.hpp"
#include "report.hpp"
#include "commit.hpp"
#include "../reads/prepare.hpp"
#include "../output/sink.hpp"

namespace nustar::star2711b {

class SingleEndRunner {
public:
    Result<RunReport> run(const SingleRunPlan&, IReadSource&, ISingleReadPreparer&, IOutputSink&);
    void cancel() noexcept;
private:
    struct Impl;
    // Owns SE pipeline workspaces, per-worker RNG, commit state, and scheduler state.
    std::shared_ptr<Impl> impl_;
};

} // namespace nustar::star2711b
