#pragma once

#include <filesystem>
#include "types.hpp"

namespace nustar::star2711b {

Result<FileIdentity> identifyFile(const std::filesystem::path& path);

} // namespace nustar::star2711b
