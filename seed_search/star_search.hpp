#pragma once

#include <cstdint>
#include <memory>
#include <vector>
#include "interface.hpp"
#include "../reference/genome.hpp"
#include "../index/suffix_array_index.hpp"

namespace nustar::star2711b {

inline constexpr std::string_view star_seed_search_domain =
    "star-2.7.11b-full-sa-sparse1-v1";

// STAR PC fields are kept as a diagnostic result so rank ranges and ordering
// can be checked without coupling the public Seed type to suffix arrays.
struct StarSeedPiece {
    std::uint64_t read_start, length, direction, multiplicity;
    std::uint64_t sa_begin, sa_end, fragment;
};
struct StarSeedCounters {
    std::uint64_t all_alignments = 0, unique_alignments = 0, multiple_alignments = 0;
    std::uint64_t min_multimaps = 0, min_multimap_length = 0, max_unique_length = 0;
    std::uint64_t no_good_piece_length = 0;
};
struct StarSeedSearchResult {
    std::vector<SplitSegment> splits;
    std::vector<StarSeedPiece> pieces;
    StarSeedCounters counters;
    std::uint64_t map_marker = 0;
};

// Immutable and safe to share between workers. Only the default STAR
// genomeSAsparseD=1 search domain is currently supported by this index.
class StarSeedSearch {
public:
    StarSeedSearch(std::shared_ptr<const SuffixArrayIndex> index,
                   std::shared_ptr<const GenomeReferenceData> reference,
                   SeedConfig policy);
    Result<StarSeedSearchResult> search(const MappingRead& read) const;
    // Reuse PC/split storage, as STAR does in its per-worker ReadAlign.
    Status search(const MappingRead& read, StarSeedSearchResult& out) const;
    const SuffixArrayIndex& index() const noexcept { return *index_; }
    const GenomeReferenceData& reference() const noexcept { return *reference_; }
private:
    struct Match {
        std::uint64_t length, begin, end, multiplicity;
    };
    Match longest(const MappingRead& read, std::uint64_t start,
                  std::uint64_t length, bool forward) const;
    template<bool WordComparison>
    std::uint64_t compare(const MappingRead& read, std::uint64_t start,
        std::uint64_t length, std::uint64_t known, std::uint64_t rank,
        bool forward, bool& query_greater) const;
    std::shared_ptr<const SuffixArrayIndex> index_;
    std::shared_ptr<const GenomeReferenceData> reference_;
    SeedConfig policy_;
    std::span<const BaseCode> genome_;
    std::uint64_t real_length_ = 0;

};

class StarSeedAlgorithmFactory final : public ISeedAlgorithmFactory {
public:
    std::string_view algorithmId() const noexcept override { return "star-2.7.11b-seed-v1"; }
    Result<std::shared_ptr<const ISeedFinder>> bind(
        SearchIndexLease, ReferenceLease,
        std::shared_ptr<const SeedSearchPlan>) const override;
};

} // namespace nustar::star2711b
