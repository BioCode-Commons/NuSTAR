#include "../runtime/two_pass.hpp"
#include "align_paired.hpp"
#include "../runtime/paired_alignment.hpp"

#include <iostream>
#include <string>

int runAlignPairedSam(const nustar::star2711b::PairedSamOptions& options,
                      bool runtime_profile) {
    using namespace nustar::star2711b;
    const auto result = options.run.twopass ? runTwoPass(options) : runPairedSam(options);
    if (const auto* failure = std::get_if<Failure>(&result)) {
        std::cerr << "alignReads: " << failure->message << '\n';
        return 1;
    }
    const auto& report = std::get<PairedSamReport>(result);
    std::cerr << "pairs=" << report.pairs << " sam_records=" << report.records << '\n';
    if (runtime_profile)
        std::cerr << "initialization_seconds total=" << report.initialization_seconds
                  << "\nmapping_wall_seconds=" << report.mapping_wall_seconds << '\n';
    return 0;
}
