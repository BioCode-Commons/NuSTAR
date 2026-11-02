#pragma once
#include "single_alignment.hpp"
#include "paired_alignment.hpp"

namespace nustar::star2711b {
Result<SingleSamReport> runTwoPass(const SingleSamOptions&);
Result<PairedSamReport> runTwoPass(const PairedSamOptions&);
}
