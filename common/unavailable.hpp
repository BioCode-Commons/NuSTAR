#pragma once

#include <string>
#include <utility>
#include "types.hpp"

namespace nustar::star2711b {

// Framework boundary: unfinished modules fail before writing compatibility output.
inline Failure unavailable(Stage stage, std::string component) {
    return Failure{ErrorKind::Config, stage, std::nullopt, 1,
                   std::move(component) + " backend is not installed"};
}

} // namespace nustar::star2711b
