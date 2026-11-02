#include "counting/feature_counter.hpp"
#include "session.hpp"

namespace nustar::star2711b::solo {
using detail::AnnotationIndex;
using detail::Context;
using detail::FeatureCounter;
Status filterCells(const Config& config, const std::filesystem::path& input,
                   const std::string& output_prefix, const std::string& log_prefix) {
    Context p;
    p.inOut->logStdOut = &std::cout;
    const auto log_parent = std::filesystem::path(log_prefix).parent_path();
    if (!log_parent.empty())
        std::filesystem::create_directories(log_parent);
    p.inOut->logMain.open(log_prefix + "Log.out");
    p.readOptions(config.arguments);
    p.runMode = "soloCellFiltering";
    p.runModeIn = {p.runMode, input.string(), output_prefix};
    p.quant = {};
    p.runDirPerm = 0700;
    p.outFileNamePrefix = log_prefix;
    p.options.initialize(&p);
    AnnotationIndex trans(p);
    FeatureCounter feature(p, nullptr, trans, -1, nullptr, nullptr);
    feature.loadRawMatrix();
    feature.cellFiltering();
    p.ownedStreams.clear();

    return std::monostate{};
}

} // namespace nustar::star2711b::solo
