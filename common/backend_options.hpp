#pragma once

#include "types.hpp"

namespace nustar::star2711b {

struct BackendOptions {
    std::string schema_id;
    std::uint32_t schema_version;
    std::vector<OptionProvenance> values; // Validated by that backend; unknown keys are errors.
};

} // namespace nustar::star2711b
