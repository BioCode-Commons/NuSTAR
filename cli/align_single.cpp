#include "../runtime/two_pass.hpp"
#include "align_single.hpp"
#include "../runtime/single_alignment.hpp"

#include <iostream>
#include <string>

int runAlignSingleSam(const nustar::star2711b::SingleSamOptions& options,
                      bool runtime_profile) {
    using namespace nustar::star2711b;
    const bool profile = options.profile_stages;
    const auto result = options.run.twopass ? runTwoPass(options) : runSingleSam(options);
    if (const auto* failure = std::get_if<Failure>(&result)) {
        std::cerr << "alignReads: " << failure->message << '\n';
        return 1;
    }
    const auto& report = std::get<SingleSamReport>(result);
    std::cerr << "reads=" << report.reads << " sam_records=" << report.records << '\n';
    if (profile) {
        std::cerr << "stages_seconds seed=" << report.seed_seconds
                  << " window=" << report.window_seconds
                  << " stitch=" << report.stitch_seconds
                  << " selection=" << report.selection_seconds
                  << " sam=" << report.sam_seconds << '\n';
    }
    if (profile || runtime_profile) {
        std::cerr << "initialization_seconds reference=" << report.reference_load_seconds
                  << " index=" << report.index_load_seconds
                  << " seed_bind=" << report.seed_bind_seconds
                  << " pipeline_setup=" << report.pipeline_setup_seconds << '\n';
        std::cerr << "mapping_wall_seconds=" << report.mapping_wall_seconds;
        if (profile)
            std::cerr << " worker_wait_seconds input=" << report.input_wait_seconds
                      << " commit=" << report.commit_wait_seconds;
        std::cerr << '\n';
    }
    if (profile) {
        std::cerr << "input_seconds read=" << report.input_read_seconds
                  << " prepare=" << report.prepare_seconds << '\n';
    }
    return 0;
}
