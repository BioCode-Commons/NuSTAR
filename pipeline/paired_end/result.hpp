#pragma once

#include "../../reads/prepare.hpp"
#include "../../stitching/candidate.hpp"
#include "overlap_state.hpp"
#include "../../selection/model.hpp"
#include "../effects.hpp"

namespace nustar::star2711b {

struct PairedReadComputation {
    PreparedReadPair read;
    CandidateSet candidates;
    OverlapState overlap;
    SelectionResult selection;
    FilterDecision filter;
    FragmentEffects effects;
};

} // namespace nustar::star2711b
