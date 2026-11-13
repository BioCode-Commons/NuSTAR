#include "../splice_junction/table_writer.hpp"
#include "alignment_scheduler.hpp"

#include "../splice_junction/star_junction.hpp"

#include <atomic>
#include <algorithm>
#include <exception>
#include <mutex>
#include <thread>

namespace nustar::star2711b {
namespace {

struct CompletedBatch {
    OutputBatch output;
    std::vector<FragmentJunctions> junctions;
    std::vector<StatsDelta> statistics;
    std::uint64_t sam_records = 0;
};

Failure runtimeFailure(std::string message) {
    return Failure{ErrorKind::Internal, Stage::Finish, {}, 1, std::move(message)};
}

struct FailureException final : std::exception {
    explicit FailureException(Failure value) : failure(std::move(value)) {}
    const char* what() const noexcept override { return failure.message.c_str(); }
    Failure failure;
};

template<class T> T unwrap(Result<T> result) {
    if (const auto* error = std::get_if<Failure>(&result))
        throw FailureException(*error);
    return std::get<T>(std::move(result));
}

void unwrap(Status result) {
    if (const auto* error = std::get_if<Failure>(&result))
        throw FailureException(*error);
}

} // namespace

void WorkerAlignmentReduction::add(BatchId id, std::vector<FragmentJunctions>& junctions,
    std::span<const StatsDelta> deltas, std::uint64_t records) {
    JunctionBatch batch{id, {}};
    std::size_t size = 0;
    for (const auto& fragment : junctions) size += fragment.ordered.size();
    batch.observations.reserve(size);
    for (const auto& fragment : junctions)
        batch.observations.insert(batch.observations.end(), fragment.ordered.begin(), fragment.ordered.end());
    for (const auto& delta : deltas) addStatistics(statistics.totals, delta);
    report.fragments += junctions.size();
    report.sam_records += records;
    batches.push_back(std::move(batch));
}

Result<ParallelAlignReport> finishConcurrentAlignment(
    std::span<WorkerAlignmentReduction> workers, IOutputSink& sink,
    const SjOutputConfig& config, const ReferenceMetadata& metadata) {
    try {
        ParallelAlignReport report;
        WorkerStatistics statistics{};
        std::vector<WorkerAlignmentReduction::JunctionBatch*> batches;
        for (auto& worker : workers) {
            report.fragments += worker.report.fragments;
            report.sam_records += worker.report.sam_records;
            addStatistics(statistics.totals, worker.statistics.totals);
            for (auto& batch : worker.batches) batches.push_back(&batch);
        }
        std::sort(batches.begin(), batches.end(), [](const auto* a, const auto* b) {
            return a->id.value < b->id.value;
        });
        JunctionAccumulator junctions(0, config);
        for (auto* batch : batches) {
            FragmentJunctions fragment{{}, std::move(batch->observations)};
            unwrap(junctions.add(fragment));
        }
        std::vector<WorkerJunctions> reduced;
        reduced.push_back(unwrap(junctions.finish()));
        auto final = unwrap(finalizeJunctions(reduced, config));
        unwrap(sink.finish(unwrap(encodeSjTable(final, metadata)), statistics));
        return report;
    } catch (const FailureException& error) { return error.failure; }
      catch (const std::exception& error) { return runtimeFailure(error.what()); }
}

Result<ParallelAlignReport> runParallelAlignment(
    IReadSource& source, IOutputSink& sink, std::span<const std::byte> header,
    const SjOutputConfig& config, const ReferenceMetadata& metadata,
    std::uint32_t threads, std::uint64_t target_bytes,
    WorkspaceFactory workspace_factory, FragmentComputer compute) {
    const bool parallel_output = sink.concurrentWrites();
    std::atomic<bool> stopped{false};
    std::mutex failure_mutex, input_mutex, output_mutex;
    std::optional<Failure> failure;
    const auto fail = [&](Failure error) {
        std::lock_guard lock(failure_mutex);
        if (!failure) failure = std::move(error);
        stopped.store(true, std::memory_order_relaxed);
    };
    std::vector<WorkerAlignmentReduction> reductions(threads);
    std::vector<std::jthread> workers;
    try {
        // Prepare the source; workers perform block reads and framing on demand.
        unwrap(source.prepareParallel(target_bytes, threads));
        unwrap(sink.begin(header));
        const bool parallel_input = source.concurrentReads();
        workers.reserve(threads);
        for (std::uint32_t id = 0; id < threads; ++id) {
            workers.emplace_back([&, id] {
                try {
                    auto workspace = unwrap(workspace_factory());
                    LegacyRandomState random{};
                    random.engine.seed(777);
                    while (!stopped.load(std::memory_order_relaxed)) {
                        std::optional<ReadBatch> task;
                        if (parallel_input) task = unwrap(source.next(target_bytes));
                        else {
                            // SAM/BAM and barcode wrappers retain their streaming contract.
                            std::lock_guard lock(input_mutex);
                            task = unwrap(source.next(target_bytes));
                        }
                        if (!task) break;
                        CompletedBatch done{OutputBatch{task->id, {}}, {}, {}, 0};
                        done.output.fragments.reserve(task->fragments.size());
                        done.junctions.reserve(task->fragments.size());
                        done.statistics.reserve(task->fragments.size());
                        for (auto& fragment : task->fragments) {
                            if (stopped.load(std::memory_order_relaxed)) return;
                            auto effects = unwrap(compute(std::move(fragment), *workspace, random));
                            done.sam_records += effects.output.records.size();
                            done.output.fragments.push_back(unwrap(sink.encode(effects.fragment, effects.output)));
                            done.junctions.push_back(std::move(effects.junctions));
                            done.statistics.push_back(effects.statistics);
                        }
                        if (parallel_output) unwrap(sink.submit(std::move(done.output)));
                        else {
                            std::lock_guard lock(output_mutex);
                            unwrap(sink.submit(std::move(done.output)));
                        }
                        reductions[id].add(task->id, done.junctions, done.statistics, done.sam_records);
                    }
                } catch (const FailureException& error) { fail(error.failure); }
                  catch (const std::exception& error) { fail(runtimeFailure(error.what())); }
            });
        }
    } catch (const FailureException& error) { fail(error.failure); }
      catch (const std::exception& error) { fail(runtimeFailure(error.what())); }
    for (auto& worker : workers) worker.join();
    if (!failure) {
        auto result = finishConcurrentAlignment(reductions, sink, config, metadata);
        if (const auto* error = std::get_if<Failure>(&result)) failure = *error;
        else return result;
    }
    source.cancel();
    sink.cancel();
    return *failure;
}


} // namespace nustar::star2711b
