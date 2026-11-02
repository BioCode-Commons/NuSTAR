#pragma once

#include <string>
#include <variant>
#include "../common/types.hpp"

namespace nustar::star2711b {

struct MappingBackendNames {
    std::string candidate_mapper, seed, window, stitch;
};
struct SingleBackendNames {
    MappingBackendNames mapping;
    std::string preparer, pipeline, selection, sam_projector;
};
struct PairedBackendNames {
    MappingBackendNames mapping, merged_single_mapping;
    std::string preparer, pipeline, overlap, selection, sam_projector;
};
struct BackendNames {
    std::string index;
    std::variant<SingleBackendNames, PairedBackendNames> layout;
};

} // namespace nustar::star2711b
