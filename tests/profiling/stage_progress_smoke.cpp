#include "../../profiling/stage_progress.hpp"

#include <fstream>
#include <filesystem>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

using namespace nustar::star2711b;

namespace {

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    for (;;) {
        const auto end = line.find('\t', begin);
        fields.push_back(line.substr(begin, end == std::string::npos ? end : end-begin));
        if (end == std::string::npos) return fields;
        begin = end+1;
    }
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    const auto stages = timedStageBit(TimedStage::Seed) |
                        timedStageBit(TimedStage::Cluster) |
                        timedStageBit(TimedStage::Stitch);
    {
        StageProgress progress(argv[1], 2, stages);
        for (int i = 0; i < 10001; ++i) {
            progress.beginRead();
            for (auto stage : {TimedStage::Seed, TimedStage::Cluster, TimedStage::Stitch}) {
                progress.beginStage(stage);
                progress.endStage(stage);
            }
            progress.endRead();
        }
        progress.finish();
    }
    std::ifstream input(argv[1]);
    if (!input) throw std::runtime_error("timing file missing");
    std::string line;
    std::vector<std::vector<std::string>> rows;
    while (std::getline(input, line)) rows.push_back(split(line));
    if (rows.size() != 3 || rows[0].size() != 9 ||
        rows[1].size() != 9 || rows[2].size() != 9 ||
        rows[1][0] != "progress" || rows[1][1] != "10000" ||
        rows[2][0] != "final" || rows[2][1] != "10001")
        throw std::runtime_error("timing progress rows differ");
    for (const auto& row : rows) {
        for (const auto& field : row)
            if (field == "NA") throw std::runtime_error("available stage written as NA");
    }

    const auto partial_path = std::string(argv[1]) + ".partial";
    {
        StageProgress partial(partial_path, 1, timedStageBit(TimedStage::Seed));
        partial.beginRead();
        partial.beginStage(TimedStage::Seed);
        partial.endStage(TimedStage::Seed);
        partial.endRead();
        partial.finish();
    }
    std::ifstream partial_input(partial_path);
    std::vector<std::vector<std::string>> partial_rows;
    while (std::getline(partial_input, line)) partial_rows.push_back(split(line));
    if (partial_rows.size() != 2 || partial_rows[1].size() != 9 ||
        partial_rows[1][0] != "final" || partial_rows[1][1] != "1" ||
        partial_rows[1][4] != "NA" || partial_rows[1][5] != "NA" ||
        partial_rows[1][6] != "NA" || partial_rows[1][8] != "NA")
        throw std::runtime_error("unavailable stages not marked NA");

    const auto disabled_path = std::string(argv[1]) + ".disabled";
    const bool existed_before = std::filesystem::exists(disabled_path);
    {
        StageProgress disabled(disabled_path, 0, stages);
        disabled.beginRead();
        disabled.beginStage(TimedStage::Seed);
        disabled.endStage(TimedStage::Seed);
        disabled.endRead();
        disabled.finish();
    }
    if (!existed_before && std::filesystem::exists(disabled_path))
        throw std::runtime_error("disabled timing created a file");

    const auto clamped_path = std::string(argv[1]) + ".clamped";
    if (setenv("NUSTAR_STAGE_TIMING_EVERY", "250", 1) != 0 ||
        unsetenv("NUSTAR_STAGE_TIMING_FILE") != 0)
        throw std::runtime_error("cannot configure timing environment");
    {
        auto clamped = StageProgress::fromEnvironment(
            clamped_path, timedStageBit(TimedStage::Seed));
        for (int i = 0; i < 10000; ++i) {
            clamped.beginRead();
            clamped.endRead();
        }
        clamped.finish();
    }
    unsetenv("NUSTAR_STAGE_TIMING_EVERY");
    std::ifstream clamped_input(clamped_path);
    std::vector<std::vector<std::string>> clamped_rows;
    while (std::getline(clamped_input, line)) clamped_rows.push_back(split(line));
    if (clamped_rows.size() != 3 || clamped_rows[1][0] != "progress" ||
        clamped_rows[1][1] != "10000" || clamped_rows[2][0] != "final" ||
        clamped_rows[2][1] != "10000")
        throw std::runtime_error("environment interval not clamped to 10000");
}
