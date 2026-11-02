#pragma once

// Reserved framework contract; see README for current implementation status.

#include <cstdint>
#include <optional>
#include <string>
#include <variant>
#include <vector>
#include "../../seed_search/delivery.hpp"
#include "../../window/model.hpp"
#include "../../pipeline/single_end/mapping_interface.hpp"
#include "../../pipeline/paired_end/overlap.hpp"
#include "../../pipeline/single_end/result.hpp"
#include "../../pipeline/paired_end/result.hpp"
#include "../../runtime/backend_names.hpp"
#include "../../index/backend.hpp"

namespace nustar::star2711b {

// Compute and trace advance only the supplied private RNG and workspace.
// Trace captures extra values but must produce identical effects and RNG state.
// Neither path commits files, SJ counts, statistics, or quality histograms.
struct StageDifference {
    FragmentId fragment;
    Stage stage;
    enum class Disposition { CompatibilityFailure, Diagnostic, NotComparable } disposition;
    std::string field_path, reference_value, candidate_value;
};
struct SnapshotComparisonPlan {
    enum class Acceptance { EquivalentSeedStreamAndFinalOutput, FinalOutput } acceptance;
    bool compare_native_trace; // Only for a comparable schema and native index representation.
};
struct ReplayIdentity {
    BaselineIdentity baseline;
    std::string input_sha256, reference_semantic_id, config_sha256;
    BackendNames resolved_backends;
    std::string main_mapping_config_sha256;
    std::optional<std::string> merged_mapping_config_sha256;
    SearchIndexIdentity search_index;
    ReferenceGeneration generation;
    std::uint32_t worker;
    BatchId batch;
    std::uint64_t worker_fragment_ordinal;
    RandomSnapshot random_before;
};
struct SeedHitSetSnapshot {
    OccurrenceHandle handle; // Rebound to a new source ID when restored; not a comparison key.
    std::vector<SeedHit> hits;
    bool complete; // False for diagnostic truncation; never sufficient for an equivalence pass.
};
struct SeedBatchSnapshot {
    ReferenceGeneration generation;
    SeedDeliveryContract delivery;
    std::vector<SplitSegment> splits;
    std::vector<Seed> seeds;
    std::vector<SeedHitSetSnapshot> occurrences; // Values only; no cursor, vtable, or index pointer.
    std::optional<BackendTrace> native_trace;
    StageMarkers markers;
};
Result<SeedBatchSnapshot> snapshotSeeds(const SeedBatch&, std::optional<std::uint64_t> hit_limit);
struct MappingStageSnapshot {
    SeedBatchSnapshot seeds;
    WindowBatch windows;
    CandidateSet candidates; // Before overlap, selection, and output projection.
};
struct MergedMappingSnapshot {
    MergedPairMapping input;
    MappingStageSnapshot mapping;
};
struct SingleFragmentSnapshot {
    MappingStageSnapshot mapping;
    SingleReadComputation result;
};
struct PairedFragmentSnapshot {
    MappingStageSnapshot mapping;
    std::optional<MergedMappingSnapshot> merged;
    PairedReadComputation result;
};
struct FragmentSnapshot {
    std::uint32_t schema_version;
    ReplayIdentity identity;
    std::variant<SingleFragmentSnapshot, PairedFragmentSnapshot> layout;
    RandomSnapshot random_after;
};
Result<std::vector<StageDifference>> compareSnapshots(
    const FragmentSnapshot& reference, const FragmentSnapshot& candidate,
    const SnapshotComparisonPlan&);

} // namespace nustar::star2711b
