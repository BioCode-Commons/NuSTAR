#pragma once

#include <vector>
#include "model.hpp"
#include "../common/compatibility_limits.hpp"

namespace nustar::star2711b {

struct ReadFilterConfig;

class ISingleReadPreparer {
public:
    virtual ~ISingleReadPreparer() = default;
    // Batch preparation preserves chunk-level clipping without a paired-mode branch.
    virtual Result<std::vector<PreparedSingleRead>> prepareSingleBatch(
        SingleReadBatch&&, const SinglePrepareConfig&, const ReadFilterConfig&,
        const CompatibilityLimits&) = 0;
};
class IPairedReadPreparer {
public:
    virtual ~IPairedReadPreparer() = default;
    virtual Result<std::vector<PreparedReadPair>> preparePairedBatch(
        PairedReadBatch&&, const PairedPrepareConfig&, const ReadFilterConfig&,
        const CompatibilityLimits&) = 0;
};

} // namespace nustar::star2711b
