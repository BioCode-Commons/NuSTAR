#pragma once

#include <variant>
#include <vector>
#include "../reads/prepare_config.hpp"
#include "../stitching/config.hpp"
#include "../pipeline/paired_end/overlap_config.hpp"
#include "../sam/config.hpp"
#include "../window/config.hpp"
#include "../stitching/score.hpp"
#include "../selection/config.hpp"
#include "../splice_junction/config.hpp"
#include "../annotation/config.hpp"
#include "../common/compatibility_limits.hpp"
#include "execution_config.hpp"

namespace nustar::star2711b {

struct SingleLayoutConfig {
    SinglePrepareConfig prepare;
    SingleStitchConfig stitch;
    SamConfig sam;
};
struct PairedLayoutConfig {
    PairedPrepareConfig prepare;
    PairedStitchConfig stitch;
    OverlapConfig overlap;
    PairedSamConfig sam;
};

struct ResolvedRunConfig {
    BaselineIdentity baseline;
    std::variant<SingleLayoutConfig, PairedLayoutConfig> layout;
    WindowConfig window;
    ScoreConfig score;
    SelectionConfig select;
    ReadFilterConfig filter;
    SjOutputConfig sj;
    AnnotationConfig annotation;
    CompatibilityLimits limits;
    RuntimeConfig runtime;
    std::vector<OptionProvenance> provenance;
    std::vector<OptionDecision> option_decisions; // Exhaustive; all Projected in a modular RunPlan.
};

} // namespace nustar::star2711b
