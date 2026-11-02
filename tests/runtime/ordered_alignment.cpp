#include "../../runtime/ordered_alignment.hpp"
#include "../../splice_junction/star_junction.hpp"

#include <atomic>
#include <chrono>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <thread>

using namespace nustar::star2711b;

namespace {
Failure injected() { return {ErrorKind::Input, Stage::Input, {}, 1, "injected"}; }
struct Workspace final : StageWorkspace {
    WorkspaceIdentity id{"test", 1, "test"};
    const WorkspaceIdentity& identity() const noexcept override { return id; }
    void reset() noexcept override {}
};
struct Source final : IReadSource {
    std::uint64_t count, current = 0;
    bool fail = false, cancelled = false;
    std::vector<std::thread::id> owners;
    std::atomic<std::uint64_t> committed{0};
    std::uint64_t limit;
    Source(std::uint64_t n, std::uint64_t bound) : count(n), owners(n), limit(bound) {}
    Result<std::optional<ReadBatch>> next(std::uint64_t) override {
        if (fail && current == 200) return injected();
        if (current == count) return std::optional<ReadBatch>{};
        if (current-committed.load() >= limit)
            throw std::runtime_error("unbounded in-flight batches");
        owners[current] = std::this_thread::get_id();
        RawFragment fragment{}; fragment.id = {0, current};
        return std::optional<ReadBatch>{ReadBatch{{current++}, {std::move(fragment)}}};
    }
    void cancel() noexcept override { cancelled = true; }
};
struct Sink final : IOutputSink {
    std::uint64_t committed = 0;
    bool fail = false, fail_finish = false, finished = false, cancelled = false;
    bool spooling=false, fail_stage=false;
    Source& source;
    explicit Sink(Source& value,bool staged) : spooling(staged),source(value) {}
    struct Chunk final : StagedOutput {
        explicit Chunk(OutputBatch value) : StagedOutput(value.batch), output(std::move(value)) { ++live; }
        ~Chunk() { --live; }
        OutputBatch output;
        static inline std::atomic<unsigned> live{0};
    };
    bool supportsStaging() const noexcept override { return spooling; }
    Result<std::unique_ptr<StagedOutput>> stage(OutputBatch&& batch) override {
        if (source.owners[batch.batch.value]!=std::this_thread::get_id())
            throw std::runtime_error("staging moved off the producing worker");
        if (fail_stage && batch.batch.value==200) return injected();
        auto chunk=std::make_unique<Chunk>(std::move(batch));
        // stage consumes the batch; callers must not use its ID afterward.
        batch.batch=BatchId{std::numeric_limits<std::uint64_t>::max()};
        return std::unique_ptr<StagedOutput>(std::move(chunk));
    }
    Status submitStaged(std::unique_ptr<StagedOutput> staged) override {
        return accept(std::move(dynamic_cast<Chunk&>(*staged).output));
    }
    Status begin(std::span<const std::byte>) override { return std::monostate{}; }
    Status submit(OutputBatch&& batch) override {
        if (source.owners[batch.batch.value] != std::this_thread::get_id())
            throw std::runtime_error("batch crossed allocation/commit workers");
        return accept(std::move(batch));
    }
    Status accept(OutputBatch&& batch) {
        if (fail && committed == 200) return injected();
        if (batch.batch.value != committed || batch.fragments.size() != 1)
            throw std::runtime_error("output lost input order or a fragment");
        ++committed;
        source.committed.store(committed);
        return std::monostate{};
    }
    Status finish(std::span<const std::byte>, const WorkerStatistics& statistics) override {
        if (fail_finish) return injected();
        if (statistics.totals.reads != committed)
            throw std::runtime_error("statistics lost a batch");
        finished = true;
        return std::monostate{};
    }
    void cancel() noexcept override { cancelled = true; }
};
void run(unsigned threads, unsigned failure, std::uint64_t count, bool staged) {
    Source source(count, 2*threads); source.fail = failure == 1;
    Sink sink(source,staged); sink.fail = failure == 2; sink.fail_finish = failure == 5;
    sink.fail_stage=failure==6;
    ReferenceMetadata metadata{{{ChromosomeId{0}, "chr", GenomeOffset{0}, 100}},
                               GenomeOffset{100}, "test"};
    auto result = runOrderedAlignment(source, sink, {}, defaultStarSjOutputConfig(),
        metadata, threads, 1,
        [failure]() -> Result<std::unique_ptr<StageWorkspace>> {
            if (failure == 4) return injected();
            return std::unique_ptr<StageWorkspace>(std::make_unique<Workspace>());
        },
        [failure, &source](RawFragment&& raw, StageWorkspace&, LegacyRandomState&) -> Result<FragmentEffects> {
            if (source.owners[raw.id.ordinal] != std::this_thread::get_id())
                throw std::runtime_error("batch crossed input/compute workers");
            if (failure == 3 && raw.id.ordinal == 200) return injected();
            // Make later batches finish first and fill the bounded in-flight
            // window, exercising bounded buffering and ordered worker commits.
            if (raw.id.ordinal%97 == 0)
                std::this_thread::sleep_for(std::chrono::milliseconds(2));
            FragmentEffects effects{}; effects.fragment = raw.id;
            effects.junctions.fragment = raw.id; effects.statistics.reads = 1;
            SamRecord record{}; record.fragment = raw.id;
            record.qname = std::to_string(raw.id.ordinal); record.sequence = "A";
            record.quality = "I"; record.rname = "*"; record.rnext = "*";
            effects.output.records.push_back(std::move(record));
            return effects;
        });
    if (Sink::Chunk::live.load()) throw std::runtime_error("staged output leaked after success or cancellation");
    if (failure) {
        if (!std::holds_alternative<Failure>(result) || !source.cancelled || !sink.cancelled)
            throw std::runtime_error("failure did not cancel all participants");
    } else if (const auto* report = std::get_if<OrderedAlignReport>(&result);
               !report || report->fragments != count || report->sam_records != count ||
               sink.committed != count || !sink.finished)
        throw std::runtime_error("ordered alignment did not finish exactly once");
}
}
int main() {
    try {
        for (const auto threads : {1U,32U}) for (const auto staged : {false,true}) {
            run(threads,0,0,staged);
            for (unsigned failure=0;failure<6;++failure) run(threads,failure,1000,staged);
            if (threads>1 && staged) run(threads,6,1000,staged);
        }
        return 0;
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
