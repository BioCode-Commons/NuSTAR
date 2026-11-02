#pragma once
#include "suffix_array_index.hpp"
#include "../reference/genome.hpp"

namespace nustar::star2711b {
// STAR sjdbBuildIndex ordering: old suffixes retain rank order, new suffixes
// are inserted at the same comparison boundaries. No whole-genome SA rebuild.
Result<SearchIndexLease> insertJunctions(const SuffixArrayIndex& old_index,
    const GenomeReferenceData& old_reference,
    std::shared_ptr<const GenomeReferenceData> reference, unsigned threads);
}
