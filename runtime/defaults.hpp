#pragma once

#include "../stitching/config.hpp"
#include "../seed_search/config.hpp"
#include "../selection/config.hpp"
#include "../sam/config.hpp"

namespace nustar::star2711b {

StitchCommonConfig defaultStarStitchCommonConfig();
SingleStitchConfig defaultStarSingleStitchConfig();
PairedStitchConfig defaultStarPairedStitchConfig();
ScoreConfig defaultStarScoreConfig();
SeedConfig defaultStarSeedConfig();
SelectionConfig defaultStarSelectionConfig();
ReadFilterConfig defaultStarReadFilterConfig();
SamConfig defaultStarSamConfig();

} // namespace nustar::star2711b
