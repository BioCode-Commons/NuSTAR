#pragma once

#include "service.hpp"

namespace nustar::star2711b {

class GtfTableWriter final : public IAnnotationTableWriter {
public:
    Status writeLegacyTables(const AnnotationCatalog&,
        std::span<const JunctionProposal>, const ReferenceMetadata&,
        const std::filesystem::path& directory) override;
};

Status writeJunctionDatabaseTables(const ReferenceView&,
    const std::filesystem::path& directory);

} // namespace nustar::star2711b
