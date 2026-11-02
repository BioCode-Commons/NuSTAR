#pragma once

#include <array>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>
#include "../common/types.hpp"

namespace nustar::star2711b {

struct InputGroup {
    std::array<std::filesystem::path, 2> mate_paths; // Only [0] is used for SE.
};
struct InputSpec {
    Layout layout;
    std::vector<InputGroup> groups;
    std::vector<std::string> readFilesCommand;
    std::optional<std::uint64_t> readMapNumber;
    bool capture_aux = false;
};

} // namespace nustar::star2711b
