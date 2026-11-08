#pragma once

#include "model.hpp"
#include "../reads/model.hpp"
#include "config.hpp"

namespace nustar::star2711b {

class IAlignmentSelector {
public:
    virtual ~IAlignmentSelector() = default;
    // Validate once when composing a pipeline, before processing any reads.
    virtual Status validateConfig(const SelectionConfig&) const { return std::monostate{}; }
    // The caller supplies a configuration accepted by validateConfig().
    virtual Status select(const MappingRead&, const CandidateSet&, const SelectionConfig&,
        LegacyRandomState&, SelectionResult& out) const = 0;
};
Result<FilterDecision> filterRead(const MappingRead&, const CandidateSet&,
    const SelectionResult&, const ReadFilterConfig&);

} // namespace nustar::star2711b
