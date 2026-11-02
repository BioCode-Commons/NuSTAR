#pragma once

#include "model.hpp"
#include "config.hpp"

namespace nustar::star2711b {

class IAnnotationLoader {
public:
    virtual ~IAnnotationLoader() = default;
    virtual Result<AnnotationCatalog> load(const AnnotationConfig&,
                                           const ReferenceView&) = 0;
    virtual Result<std::vector<JunctionProposal>> deriveJunctions(
        const AnnotationCatalog&, const ReferenceView&) = 0;
};
class IAnnotationTableWriter {
public:
    virtual ~IAnnotationTableWriter() = default;
    virtual Status writeLegacyTables(const AnnotationCatalog&,
        std::span<const JunctionProposal>, const ReferenceMetadata&,
        const std::filesystem::path& directory) = 0;
};
class IJunctionPreparer {
public:
    virtual ~IJunctionPreparer() = default;
    virtual Result<PreparedJunctionDatabase> prepare(
        std::span<const JunctionProposal>, const ReferenceView&,
        const AnnotationConfig&) = 0;
};

} // namespace nustar::star2711b
