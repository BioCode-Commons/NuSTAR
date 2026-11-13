#pragma once

#include "../pipeline/effects.hpp"

#include "../reads/raw.hpp"
#include "../output/sink.hpp"
#include "../splice_junction/interface.hpp"
#include "../common/workspace.hpp"

#include <functional>

namespace nustar::star2711b {

struct ParallelAlignReport {
    std::uint64_t fragments = 0, sam_records = 0;
};

// Worker-owned reduction state for concurrent file output. Preserve SJ
// observation order at finalization without serializing the mapping/I/O loop.
struct WorkerAlignmentReduction {
    struct JunctionBatch {
        BatchId id;
        std::vector<JunctionObservation> observations;
    };
    WorkerStatistics statistics{};
    ParallelAlignReport report;
    std::vector<JunctionBatch> batches;
    void add(BatchId, std::vector<FragmentJunctions>&, std::span<const StatsDelta>,
             std::uint64_t sam_records);
};
Result<ParallelAlignReport> finishConcurrentAlignment(
    std::span<WorkerAlignmentReduction>, IOutputSink&, const SjOutputConfig&,
    const ReferenceMetadata&);

using WorkspaceFactory = std::function<Result<std::unique_ptr<StageWorkspace>>() >;
using FragmentComputer = std::function<Result<FragmentEffects>(
    RawFragment&&, StageWorkspace&, LegacyRandomState&)>;

// Workers submit completed batches immediately. Non-concurrent sinks use a
// writer mutex, without waiting for any input batch number. SJ observations
// are reduced only after workers join; they do not gate output or input.
Result<ParallelAlignReport> runParallelAlignment(
    IReadSource&, IOutputSink&, std::span<const std::byte> sam_header,
    const SjOutputConfig&, const ReferenceMetadata&, std::uint32_t threads,
    std::uint64_t target_batch_bytes, WorkspaceFactory, FragmentComputer);

} // namespace nustar::star2711b
