#include "paired_runner.hpp"
#include "../common/unavailable.hpp"

namespace nustar::star2711b {

Result<RunReport> PairedEndRunner::run(
    const PairedRunPlan&, IReadSource&, IPairedReadPreparer&, IOutputSink&) {
    return unavailable(Stage::Plan, "paired-read runner");
}

void PairedEndRunner::cancel() noexcept {}

} // namespace nustar::star2711b
