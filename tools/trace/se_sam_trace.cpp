#include "../../runtime/single_alignment.hpp"

#include <iostream>
#include <string>
#include <string_view>

int main(int argc, char** argv) {
    if ((argc != 4 && argc != 5) ||
        (argc == 5 && std::string_view(argv[4]) != "--profile-stages")) {
        std::cerr << "Usage: se-sam-trace NUSTAR_INDEX_DIR READS.fastq OUTPUT.sam "
                     "[--profile-stages]\n";
        return 2;
    }
    using namespace nustar::star2711b;
    const std::string command = std::string("se-sam-trace ") + argv[1] +
        " " + argv[2] + " " + argv[3];
    const auto result = runSingleSam(SingleSamOptions{
        argv[1], argv[2], argv[3], command, argc == 5});
    if (const auto* failure = std::get_if<Failure>(&result)) {
        std::cerr << "se-sam-trace: " << failure->message << '\n';
        return 1;
    }
    const auto& report = std::get<SingleSamReport>(result);
    std::cerr << "se-sam-trace reads=" << report.reads
              << " records=" << report.records << '\n';
    if (argc == 5)
        std::cerr << "stages_seconds seed=" << report.seed_seconds
                  << " window=" << report.window_seconds
                  << " stitch=" << report.stitch_seconds
                  << " selection=" << report.selection_seconds
                  << " sam=" << report.sam_seconds << '\n';
    return 0;
}
