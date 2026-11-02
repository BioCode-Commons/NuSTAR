#pragma once

#include <array>
#include <cstdint>
#include <optional>
#include <vector>
#include "../common/types.hpp"
#include "../common/stage.hpp"

namespace nustar::star2711b {

enum class BoundaryKind { Intron, Deletion, Insertion, MateGap };
struct AlignmentBlock {
    ReadOffset read_start;
    GenomeOffset genome_start;
    std::uint64_t length;
    MateId mate;
    std::optional<JunctionId> sjdb; // EX_sjA preserves the seed's synthetic SJ identity.
};
struct AlignmentBoundary {
    BoundaryKind kind;
    Motif motif; // Meaningful for Intron only; noncanonical remains an intron.
    std::uint8_t annotation; // Preserve the STAR value, not a recomputed GTF lookup.
    Strand strand;
    std::uint64_t shift_left, shift_right;
};
struct AlignmentMetrics {
    Score score;
    std::uint64_t matches, mismatches, mapped_length, extension_length;
    std::uint64_t gap_count, gap_bases, deletion_count, deletion_bases;
    std::uint64_t insertion_count, insertion_bases, unique_seeds, anchor_seeds;
    std::array<std::uint64_t, 3> intron_motif_counts;
    Strand motif_strand;
    bool has_splice;
};
struct AlignmentCandidate {
    CandidateId id;
    WindowId window;
    ChromosomeId chromosome;
    Orientation orientation, original_read_orientation;
    std::optional<MateId> single_mate; // Null means that both mates contribute.
    ReadOffset read_start, original_read_start;
    std::uint64_t read_length; // Legacy rLength; do not recompute from span or mapped_length.
    GenomeOffset genome_start;
    std::uint64_t genome_length, chromosome_start;
    std::vector<AlignmentBlock> blocks;
    std::vector<AlignmentBoundary> boundaries; // Exactly blocks.size() - 1 for nonempty candidates.
    AlignmentMetrics metrics;
};
struct WindowCandidates {
    WindowId window;
    std::vector<CandidateId> ordered;
};
struct CandidateSet {
    ReferenceGeneration generation;
    std::vector<AlignmentCandidate> candidates; // ID-indexed storage; no pointer-based identity.
    std::vector<WindowCandidates> windows; // Preserve trAll/nWinTr visitation order.
    std::optional<CandidateId> best; // Legacy trBest identity, even after later reordering.
    StageMarkers markers;
};

} // namespace nustar::star2711b
