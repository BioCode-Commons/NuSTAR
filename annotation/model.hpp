#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include "../common/types.hpp"
#include "../reference/view.hpp"

namespace nustar::star2711b {

// GTF models preserve encounter IDs and exon rows; they are not alignment candidates.
struct GeneModel {
    GeneId id; // Zero-based first-encounter ID; +1 only at the legacy GTF-SJ file boundary.
    std::string external_id, name, biotype;
};
struct ExonModel {
    TranscriptId transcript;
    GeneId gene;
    ChromosomeId chromosome;
    GenomeRange interval;
    std::uint64_t source_line, encounter_order;
};
struct TranscriptModel {
    TranscriptId id;
    std::string external_id;
    Strand strand; // The first encountered transcript strand, as in STAR.
    std::vector<std::uint64_t> exon_indices; // No implicit union of overlapping exons.
};
struct AnnotationCatalog {
    std::vector<GeneModel> genes;
    std::vector<TranscriptModel> transcripts;
    std::vector<ExonModel> exons;
    std::vector<Diagnostic> diagnostics;
    bool active; // False when the GTF path is absent or effective overhang is zero.
};
enum class JunctionOrigin { ExistingIndex, JunctionFile, Gtf };
struct JunctionProposal {
    ChromosomeId chromosome;
    GenomeRange intron;
    Strand strand;
    JunctionOrigin origin;
    std::uint8_t legacy_priority; // Existing index 30; GTF 20; junction file 10.
    std::uint64_t source_order;
    std::vector<GeneId> genes; // Sorted unique IDs for GTF-derived junction metadata.
};
struct PreparedJunctionDatabase {
    std::vector<SjdbEntry> entries;
    std::vector<std::uint64_t> winning_proposal_indices;
    std::uint64_t overhang;
    std::vector<Diagnostic> diagnostics;
};

} // namespace nustar::star2711b
