#pragma once

#include "model.hpp"
#include "../reads/model.hpp"
#include "config.hpp"

namespace nustar::star2711b {

class IAlignmentSelector {
public:
    virtual ~IAlignmentSelector() = default;
    virtual Status select(const MappingRead&, const CandidateSet&, const SelectionConfig&,
        LegacyRandomState&, SelectionResult& out) const = 0;
};
Result<FilterDecision> filterRead(const MappingRead&, const CandidateSet&,
    const SelectionResult&, const ReadFilterConfig&);

} // namespace nustar::star2711b
