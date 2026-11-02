#pragma once

#include <utility>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <variant>
#include <vector>
#include "../common/stage.hpp"
#include "../reference/view.hpp"
#include "../reads/model.hpp"

namespace nustar::star2711b {

enum class HitOrder { Star2711b, BackendNative };
enum class SeedContract { Star2711bEquivalent, AlternativeStrategy };

// Delivery semantics are independent of algorithm tuning and native query capabilities.
struct SeedDeliveryContract {
    SeedContract strategy;
    HitOrder hit_order;
    bool exact_multiplicity;
    bool complete_occurrences; // Complete for each retained seed within the declared search domain.
    bool known_junction_projection;
    std::uint32_t hit_schema_version;
    std::string search_domain_id; // Versioned semantic fingerprint, not an index artifact hash.
};

// Public snapshots own values; backend scratch may use packed arrays and arenas.
struct SplitSegment { ReadRange read; MateId mate; };
enum class CountPrecision { Exact, LowerBound };
struct MatchMultiplicity {
    std::uint64_t value;
    CountPrecision precision; // A capped lower bound cannot be used as an exact repeat count.
};
struct OccurrenceHandle {
    OccurrenceSourceId source;
    OccurrenceSetId set;
};
struct OccurrenceSourceIdentity {
    OccurrenceSourceId source;
    FragmentId fragment;
    ReferenceGeneration generation;
    std::string reference_semantic_id;
    SeedDeliveryContract delivery;
};
struct Seed {
    ReadRange read; // Forward combined-read coordinates, independent of search direction.
    MateId mate;
    MatchMultiplicity multiplicity; // Legacy-compatible search-domain count, not projected hit count.
    OccurrenceHandle occurrences; // A source-scoped handle, never a suffix-array interval.
};
struct SeedHitBlock {
    ReadOffset read_start; // In the hit's oriented combined-read coordinate system.
    GenomeOffset genome_start;
    std::uint64_t length;
    MateId mate;
    std::optional<JunctionId> junction;
};
struct SeedHit {
    ChromosomeId chromosome;
    Orientation orientation;
    std::uint64_t visit_ordinal; // Zero-based projected stream position, never a physical SA rank.
    std::vector<SeedHitBlock> blocks; // One genomic block or ordered blocks across a known junction.
};
struct SeedHitView {
    ChromosomeId chromosome;
    Orientation orientation;
    std::uint64_t visit_ordinal;
    std::span<const SeedHitBlock> blocks; // Valid only during the callback.
};
struct SeedHitPage {
    std::vector<SeedHit> hits;
    bool end;
};
class ISeedHitCursor {
public:
    virtual ~ISeedHitCursor() = default;
    // Each cursor retains its own source/index leases and mutable iteration position.
    // Produce at most max_hits logical hits; empty pages are allowed only at EOF.
    virtual Status next(std::uint64_t max_hits, SeedHitPage& out) = 0; // Reuse caller-owned capacity.
};
class ISeedOccurrenceSource {
public:
    virtual ~ISeedOccurrenceSource() = default;
    virtual const OccurrenceSourceIdentity& identity() const noexcept = 0;
    // Fresh cursors support the anchor pass, distribution pass, and shadow replay.
    virtual Result<std::unique_ptr<ISeedHitCursor>> open(
        OccurrenceHandle) const = 0; // Caller supplies this source's handle; order is fixed by delivery.
    // The source owns every index/query buffer needed by its cursors; no worker scratch aliases.
    // The default implementation adapts cursors. Hot backends can override it
    // to deliver stack-backed blocks without changing the public hit format.
    virtual Status visit(OccurrenceHandle handle,
        const std::function<bool(const SeedHitView&)>& callback) const {
        auto opened = open(handle);
        if (std::holds_alternative<Failure>(opened))
            return std::get<Failure>(std::move(opened));
        auto cursor = std::get<std::unique_ptr<ISeedHitCursor>>(std::move(opened));
        SeedHitPage page{};
        do {
            auto status = cursor->next(128, page);
            if (std::holds_alternative<Failure>(status)) return status;
            for (const auto& hit : page.hits)
                if (!callback(SeedHitView{hit.chromosome, hit.orientation,
                                          hit.visit_ordinal, hit.blocks}))
                    return std::monostate{};
        } while (!page.end);
        return std::monostate{};
    }
    // Ordered oriented blocks, without chromosome/ordinal/hit-span packaging.
    // Synchronous: each block reference is valid only during its callback.
    // Stop also applies between donor and acceptor of the same occurrence.
    virtual Status visitBlocks(OccurrenceHandle handle,
        const std::function<bool(const SeedHitBlock&, Orientation)>& callback) const {
        return visit(handle, [&callback](const SeedHitView& hit) {
            for (const auto& block : hit.blocks)
                if (!callback(block, hit.orientation)) return false;
            return true;
        });
    }
    // Anchor pass needs only oriented genomic starts. A junction contributes
    // donor then acceptor; false stops the entire occurrence stream immediately.
    // Backends may avoid preparing read coordinates and block metadata here.
    virtual Status visitAnchors(OccurrenceHandle handle,
        const std::function<bool(GenomeOffset, Orientation)>& callback) const {
        return visit(handle, [&callback](const SeedHitView& hit) {
            for (const auto& block : hit.blocks)
                if (!callback(block.genome_start, hit.orientation)) return false;
            return true;
        });
    }

    // Synchronous borrowed batches amortize type erasure over several hits.
    // A false callback discards the remaining stream. An override that prepares
    // hits ahead of consumption must preserve its per-hit visitor's early-stop
    // and error ordering; the default adapter inherits visit()'s behavior.
    virtual Status visitBatch(OccurrenceHandle handle,
        const std::function<bool(std::span<const SeedHitView>)>& callback) const {
        return visit(handle,[&callback](const SeedHitView& hit) {
            return callback(std::span(&hit,1));
        });
    }

};
struct BackendTrace {
    std::string schema_id;
    std::uint32_t schema_version;
    std::vector<std::byte> payload; // Optional pointer-free diagnostic data, not a core input.
};
struct SeedBatch {
    ReferenceGeneration generation;
    std::vector<SplitSegment> splits;
    SeedDeliveryContract delivery;
    std::vector<Seed> seeds; // Ordered logical seeds; no storage-format fields.
    std::shared_ptr<const ISeedOccurrenceSource> occurrences;
    std::optional<BackendTrace> native_trace; // STAR PC/summary data stays in a versioned trace.
    StageMarkers markers;
};

} // namespace nustar::star2711b
