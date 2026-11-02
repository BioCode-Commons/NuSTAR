#pragma once

#include "service.hpp"

namespace nustar::star2711b {

class GenomeJunctionPreparer final : public IJunctionPreparer {
public:
    Result<PreparedJunctionDatabase> prepare(
        std::span<const JunctionProposal>, const ReferenceView&,
        const AnnotationConfig&) override;
};

// Reserved STAR compatibility contract; no implementation is installed.
namespace star_compat {
// Reserved adapter declaration. The product path uses reference/synthetic_junction.hpp.
Result<std::vector<BaseCode>> makeSyntheticSequence(
    const PreparedJunctionDatabase&, const ReferenceView&);
} // namespace star_compat

} // namespace nustar::star2711b
