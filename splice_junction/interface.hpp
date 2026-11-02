#pragma once

#include <cstdint>
#include <memory>
#include <span>
#include <vector>
#include "../stitching/candidate.hpp"
#include "../selection/model.hpp"
#include "../reads/model.hpp"
#include "config.hpp"

namespace nustar::star2711b {

// Junction reduction intentionally uses (global start, gap), not strand or motif, as its key.
struct JunctionKey { GenomeOffset start; std::uint32_t gap; };
struct JunctionObservation {
    JunctionKey key;
    Strand strand;
    Motif motif;
    std::uint8_t annotation;
    std::uint16_t overhang; // STAR narrows min(adjacent block lengths) before reduction.
    bool unique_fragment;
};
struct FragmentJunctions {
    FragmentId fragment;
    std::vector<JunctionObservation> ordered; // First-seen key order; max overhang on duplicates.
};
struct JunctionAggregate {
    JunctionKey key;
    Strand strand;
    Motif motif;
    std::uint8_t annotation;
    std::uint32_t unique_count, multiple_count;
    std::uint16_t max_left_overhang, max_right_overhang;
};
struct WorkerJunctions {
    std::uint32_t worker;
    std::vector<JunctionAggregate> sorted; // Sorted by key after legacy-equivalent collapse.
};
struct JunctionFilterTrace {
    JunctionKey key;
    bool support_pass, donor_distance_pass, acceptor_distance_pass, accepted;
};
struct FinalJunctions {
    std::vector<JunctionAggregate> accepted; // Final start/gap order; no absent GTF-only rows.
    std::vector<JunctionFilterTrace> trace; // Optional diagnostics, never written to SJ.out.tab.
};
Result<FragmentJunctions> collectJunctions(const MappingRead&, const CandidateSet&,
    const SelectionResult&, const FilterDecision&, const SjOutputConfig&);
class JunctionAccumulator {
public:
    JunctionAccumulator(std::uint32_t worker, const SjOutputConfig&);
    ~JunctionAccumulator();
    Status add(const FragmentJunctions&); // Called exactly once per committed fragment.
    Result<WorkerJunctions> finish(); // Consumes the accumulator; no subsequent add.
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
Result<FinalJunctions> finalizeJunctions(std::span<const WorkerJunctions>,
    const SjOutputConfig&); // Support filtering precedes donor/acceptor distance filtering.

} // namespace nustar::star2711b
