#pragma once

// Reserved framework contract; see README for current implementation status.

#include <memory>
#include "paired_plan.hpp"
#include "report.hpp"
#include "commit.hpp"
#include "../reads/prepare.hpp"
#include "../output/sink.hpp"

namespace nustar::star2711b {

class PairedEndRunner {
public:
    Result<RunReport> run(const PairedRunPlan&, IReadSource&, IPairedReadPreparer&, IOutputSink&);
    void cancel() noexcept;
private:
    struct Impl;
    // Owns PE pipeline workspaces, private overlap scratch, RNG, and commit/scheduler state.
    std::shared_ptr<Impl> impl_;
};

} // namespace nustar::star2711b
