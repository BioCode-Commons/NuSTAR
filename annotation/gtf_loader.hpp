#pragma once

#include "service.hpp"

namespace nustar::star2711b {

class GtfAnnotationLoader final : public IAnnotationLoader {
public:
    Result<AnnotationCatalog> load(const AnnotationConfig&, const ReferenceView&) override;
    Result<std::vector<JunctionProposal>> deriveJunctions(
        const AnnotationCatalog&, const ReferenceView&) override;
};

} // namespace nustar::star2711b
