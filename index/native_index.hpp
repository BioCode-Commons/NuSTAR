#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include "backend.hpp"

namespace nustar::star2711b {

// This namespace is an adapter detail. Core/window/stitch never include its native types.
namespace star_compat {
using SaRank = U64<struct SaRankTag>;
struct SaRange { SaRank begin, end; }; // Half-open; convert the legacy inclusive endpoint.
struct PackedArrayView {
    std::span<const std::byte> storage;
    std::uint64_t element_count;
    std::uint32_t bits_per_element;
    std::uint64_t value(std::uint64_t index) const;
};
struct PrefixIndexView {
    PackedArrayView entries;
    std::span<const std::uint64_t> level_starts;
    std::uint32_t prefix_bases;
    std::uint64_t absent_mask, n_mask;
};
struct SyntheticJunctionLayout {
    GenomeOffset synthetic_begin;
    std::uint64_t stride;
    std::span<const JunctionId> physical_to_logical;
};
struct StarSaView {
    std::span<const BaseCode> physical_genome;
    PackedArrayView suffix_array;
    PrefixIndexView prefix_index;
    SyntheticJunctionLayout synthetic_junctions;
    std::uint64_t genome_length, sa_sparse_step;
    std::uint32_t strand_bit;
    std::uint64_t strand_mask;
};
class IStarSaIndex : public ISearchIndex {
public:
    virtual StarSaView nativeView() const = 0;
};
} // namespace star_compat

} // namespace nustar::star2711b
