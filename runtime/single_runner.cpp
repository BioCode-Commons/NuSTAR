#include "single_runner.hpp"
#include "../common/unavailable.hpp"

namespace nustar::star2711b {

Result<RunReport> SingleEndRunner::run(
    const SingleRunPlan&, IReadSource&, ISingleReadPreparer&, IOutputSink&) {
    return unavailable(Stage::Plan, "single-read runner");
}

void SingleEndRunner::cancel() noexcept {}

} // namespace nustar::star2711b
