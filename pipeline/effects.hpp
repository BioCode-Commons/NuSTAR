#pragma once

#include "../splice_junction/interface.hpp"
#include "../sam/model.hpp"
#include "../statistics/interface.hpp"

namespace nustar::star2711b {

struct FragmentEffects {
    FragmentId fragment;
    FragmentJunctions junctions;
    ProjectedOutput output;
    StatsDelta statistics;
    QualityHistogramDelta quality; // Dense counters live only in WorkerStatistics.
};

} // namespace nustar::star2711b
