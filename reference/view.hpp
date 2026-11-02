#pragma once

#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <vector>
#include "../common/types.hpp"

namespace nustar::star2711b {

// Reference metadata outlives all aligners, SAM headers, and SJ writers.
struct Chromosome {
    ChromosomeId id;
    std::string name;
    GenomeOffset start;
    std::uint64_t length;
};
// These offsets describe a stable logical coordinate system, not physical storage.
struct ReferenceMetadata {
    std::vector<Chromosome> chromosomes; // Preserve the compatibility reference order.
    GenomeOffset real_genome_end; // Logical end including compatibility padding.
    std::string coordinate_schema; // Changes require explicit coordinate translation.
};
struct SjdbEntry {
    JunctionId id;
    GenomeRange lookup_intron;
    GenomeRange reported_intron;
    Strand strand;
    Motif motif;
    std::uint8_t shift_left, shift_right;
    GenomeOffset donor_sequence_begin, acceptor_sequence_begin;
};
struct JunctionDatabaseView {
    std::span<const SjdbEntry> entries; // Semantic IDs, independent of physical index order.
    std::uint64_t overhang; // Defines the compatibility search domain, not its storage layout.
};
class ISequenceLease {
public:
    virtual ~ISequenceLease() = default;
};
struct SequenceWindow {
    GenomeRange range;
    std::span<const BaseCode> bases;
    std::shared_ptr<const ISequenceLease> lease; // Keeps mapped or decoded bytes alive.
};
class IReferenceSequence {
public:
    virtual ~IReferenceSequence() = default;
    // Caller supplies an in-range contiguous region; inner loops use its span directly.
    virtual Result<SequenceWindow> pin(GenomeRange) const = 0;
};
struct ReferenceView {
    ReferenceGeneration generation;
    std::string_view semantic_id; // Sequence, coordinates, junction catalog, and annotation policy.
    const ReferenceMetadata* metadata;
    const IReferenceSequence* sequence;
    JunctionDatabaseView junctions;
};
class IReferenceData {
public:
    virtual ~IReferenceData() = default;
    virtual ReferenceView view() const = 0;
};
using ReferenceLease = std::shared_ptr<const IReferenceData>;

} // namespace nustar::star2711b
