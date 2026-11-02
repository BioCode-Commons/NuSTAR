#pragma once

#include <cstdint>
#include <vector>
#include "../common/types.hpp"
#include "backend_names.hpp"
#include "../statistics/interface.hpp"

namespace nustar::star2711b {

struct RunReport {
    WorkerStatistics statistics;
    std::uint64_t accepted_junctions;
    std::vector<FileIdentity> outputs;
    BackendNames resolved_backends;
};

} // namespace nustar::star2711b
