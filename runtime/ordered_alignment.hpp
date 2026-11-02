#pragma once

#include "../pipeline/effects.hpp"

#include "../reads/raw.hpp"
#include "../output/sink.hpp"
#include "../splice_junction/interface.hpp"
#include "../common/workspace.hpp"

#include <functional>

namespace nustar::star2711b {

struct OrderedAlignReport {
    std::uint64_t fragments = 0, sam_records = 0;
};

using WorkspaceFactory = std::function<Result<std::unique_ptr<StageWorkspace>>() >;
using FragmentComputer = std::function<Result<FragmentEffects>(
    RawFragment&&, StageWorkspace&, LegacyRandomState&)>;

// Workers own input, computation and SAM-buffer lifetimes. File sinks can
// spool chunks concurrently and concatenate them in input order, with at most
// 2*threads batches in flight. Other sinks retain ordered worker commits.
// Auxiliary statistics and junctions always reduce in exact input order.
Result<OrderedAlignReport> runOrderedAlignment(
    IReadSource&, IOutputSink&, std::span<const std::byte> sam_header,
    const SjOutputConfig&, const ReferenceMetadata&, std::uint32_t threads,
    std::uint64_t target_batch_bytes, WorkspaceFactory, FragmentComputer);

} // namespace nustar::star2711b
