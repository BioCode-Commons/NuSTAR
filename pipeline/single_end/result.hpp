#pragma once

#include "../../reads/prepare.hpp"
#include "../../stitching/candidate.hpp"
#include "../../selection/model.hpp"
#include "../effects.hpp"

namespace nustar::star2711b {

struct SingleReadComputation {
    PreparedSingleRead read;
    CandidateSet candidates;
    SelectionResult selection;
    FilterDecision filter;
    FragmentEffects effects;
};

} // namespace nustar::star2711b
