#pragma once

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>
#include "../common/types.hpp"

namespace nustar::star2711b {

struct AnnotationConfig {
    std::optional<std::filesystem::path> sjdbGTFfile;
    std::string sjdbGTFchrPrefix, sjdbGTFfeatureExon;
    std::string sjdbGTFtagExonParentTranscript, sjdbGTFtagExonParentGene;
    std::vector<std::string> sjdbGTFtagExonParentGeneName, sjdbGTFtagExonParentGeneType;
    std::vector<std::filesystem::path> sjdbFileChrStartEnd;
    std::uint64_t sjdbOverhang, limitSjdbInsertNsj;
    std::string sjdbInsertSave;
};

} // namespace nustar::star2711b
