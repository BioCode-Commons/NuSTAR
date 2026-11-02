#pragma once

// Reserved framework contract; see README for current implementation status.

#include <cstdint>
#include <memory>
#include "../pipeline/effects.hpp"
#include "../output/sink.hpp"
#include "../splice_junction/interface.hpp"
#include "../statistics/interface.hpp"

namespace nustar::star2711b {

class WorkerCommitState {
public:
    explicit WorkerCommitState(std::uint32_t worker);
    ~WorkerCommitState();
    WorkerCommitState(const WorkerCommitState&) = delete;
    WorkerCommitState& operator=(const WorkerCommitState&) = delete;
    WorkerJunctions finishJunctions();
    WorkerStatistics statistics() const;
private:
    struct Impl;
    std::unique_ptr<Impl> impl_; // SJ, statistics, and commit progress only; no mapping algorithms.
    friend Status commitFragment(FragmentEffects&&, WorkerCommitState&, OutputBatch&);
};
// Move-only submission is enforced by worker state; repeated fragment commits are rejected.
Status commitFragment(FragmentEffects&&, WorkerCommitState&, OutputBatch&);

} // namespace nustar::star2711b
