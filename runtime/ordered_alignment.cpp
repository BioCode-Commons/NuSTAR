#include "../splice_junction/table_writer.hpp"
#include "ordered_alignment.hpp"

#include "../splice_junction/star_junction.hpp"

#include <atomic>
#include <condition_variable>
#include <exception>
#include <mutex>
#include <map>
#include <stdexcept>
#include <thread>

namespace nustar::star2711b {
namespace {

struct CompletedBatch {
    OutputBatch output;
    std::vector<FragmentJunctions> junctions;
    std::vector<StatsDelta> statistics;
    std::uint64_t sam_records = 0;
    std::unique_ptr<StagedOutput> staged;
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

struct Shared {
    std::mutex input_mutex, commit_mutex;
    std::condition_variable committed;
    std::atomic<bool> stopped{false};
    std::optional<Failure> error; // Protected by commit_mutex until workers join.
    std::uint64_t published = 0;
    std::atomic<std::uint64_t> next_commit{0};
    bool draining=false; // One worker owns sink/auxiliary reduction at a time.
    std::map<std::uint64_t,CompletedBatch> ready; // Protected by commit_mutex.
    bool input_done = false; // Protected by input_mutex.

    void fail(Failure failure) {
        {
            std::lock_guard lock(commit_mutex);
            if (!error) error = std::move(failure);
            stopped.store(true, std::memory_order_relaxed);
        }
        committed.notify_all();
    }
};

} // namespace

Result<OrderedAlignReport> runOrderedAlignment(
    IReadSource& source, IOutputSink& sink, std::span<const std::byte> sam_header,
    const SjOutputConfig& sj_config, const ReferenceMetadata& metadata,
    std::uint32_t threads, std::uint64_t target_batch_bytes,
    WorkspaceFactory workspace_factory, FragmentComputer compute) {
    if (threads == 0 || threads > 256 || target_batch_bytes == 0)
        return Failure{ErrorKind::Config, Stage::Plan, {}, 1,
                       "invalid worker count or FASTQ batch size"};
    if (const auto status = sink.begin(sam_header);
        std::holds_alternative<Failure>(status))
        return std::get<Failure>(status);

    const bool spooling=threads>1 && sink.supportsStaging();
    Shared state;
    OrderedAlignReport report;
    WorkerStatistics statistics{};
    std::unique_ptr<JunctionAccumulator> junctions;
    std::vector<std::jthread> workers;
    try {
        junctions = std::make_unique<JunctionAccumulator>(0, sj_config);
        workers.reserve(threads);
        for (std::uint32_t worker = 0; worker < threads; ++worker) {
            workers.emplace_back([&] {
                try {
                    auto workspace = unwrap(workspace_factory());
                    LegacyRandomState random{};
                    random.engine.seed(777);
                    for (;;) {
                        // Allocate, consume and release a batch in the same
                        // worker. Transferring millions of FASTQ strings from
                        // one producer to many freeing workers contends on the
                        // producer's malloc arena, even with a lock-free queue.
                        std::optional<ReadBatch> task;
                        {
                            std::lock_guard lock(state.input_mutex);
                            if (state.stopped.load(std::memory_order_relaxed) ||
                                state.input_done) break;
                            if (spooling && state.published-state.next_commit.load(std::memory_order_relaxed)>=2*threads) {
                                std::unique_lock commit_lock(state.commit_mutex);
                                state.committed.wait(commit_lock,[&] {
                                    return state.stopped.load(std::memory_order_relaxed) ||
                                           state.published-state.next_commit<2*threads;
                                });
                                if (state.stopped.load(std::memory_order_relaxed)) break;
                            }
                            task = unwrap(source.next(target_batch_bytes));
                            if (!task) {
                                state.input_done = true;
                                break;
                            }
                            if (task->id.value != state.published)
                                throw std::runtime_error("missing completed batch");
                            ++state.published;
                        }
                        CompletedBatch done{OutputBatch{task->id, {}}, {}, {}, 0};
                        done.output.fragments.reserve(task->fragments.size());
                        done.junctions.reserve(task->fragments.size());
                        done.statistics.reserve(task->fragments.size());
                        for (auto& fragment : task->fragments) {
                            if (state.stopped.load(std::memory_order_relaxed)) return;
                            auto effects = unwrap(compute(std::move(fragment),
                                                          *workspace, random));
                            done.sam_records += effects.output.records.size();
                            done.output.fragments.push_back(unwrap(sink.encode(
                                effects.fragment, effects.output)));
                            done.junctions.push_back(std::move(effects.junctions));
                            done.statistics.push_back(effects.statistics);
                        }
                        if (spooling) done.staged=unwrap(sink.stage(std::move(done.output)));
                        {
                            std::unique_lock lock(state.commit_mutex);
                            auto apply=[&](CompletedBatch& batch) {
                                for (std::size_t i=0;i<batch.junctions.size();++i) {
                                    unwrap(junctions->add(batch.junctions[i]));
                                    addStatistics(statistics.totals,batch.statistics[i]);
                                }
                                report.fragments+=batch.junctions.size();
                                report.sam_records+=batch.sam_records;
                                if (spooling) unwrap(sink.submitStaged(std::move(batch.staged)));
                                else unwrap(sink.submit(std::move(batch.output)));
                                ++state.next_commit;
                            };
                            if (spooling) {
                                if (state.stopped.load(std::memory_order_relaxed)) return;
                                state.ready.emplace(task->id.value,std::move(done));
                                if (!state.draining) {
                                    state.draining=true;
                                    while (!state.stopped.load(std::memory_order_relaxed)) {
                                        auto next=state.ready.find(state.next_commit.load());
                                        if (next==state.ready.end()) break;
                                        auto batch=std::move(next->second);
                                        state.ready.erase(next);
                                        // File concatenation must not block other workers
                                        // from publishing completed chunks or reading input.
                                        lock.unlock(); apply(batch); lock.lock();
                                        state.committed.notify_all();
                                    }
                                    state.draining=false;
                                }
                            } else {
                                state.committed.wait(lock,[&] {
                                    return state.stopped.load(std::memory_order_relaxed) ||
                                           done.output.batch.value==state.next_commit;
                                });
                                if (state.stopped.load(std::memory_order_relaxed)) return;
                                apply(done);
                            }
                        }
                        // Waiters have different batch predicates: notify_one
                        // could wake the wrong batch and leave the next asleep.
                        state.committed.notify_all();
                        // FASTQ/SAM allocations die on their producing worker. Only
                        // bounded chunk handles and auxiliary deltas can be queued.
                    }
                } catch (const FailureException& error) {
                    state.fail(error.failure);
                } catch (const std::exception& error) {
                    state.fail(runtimeFailure(error.what()));
                }
            });
        }
    } catch (const std::exception& error) {
        state.fail(runtimeFailure(error.what()));
    }
    for (auto& worker : workers) worker.join();
    if (!state.error) {
        try {
            if (state.next_commit != state.published)
                throw std::runtime_error("missing completed batch");
            std::vector<WorkerJunctions> reduced;
            reduced.push_back(unwrap(junctions->finish()));
            const auto final = unwrap(finalizeJunctions(reduced, sj_config));
            const auto bytes = unwrap(encodeSjTable(final, metadata));
            unwrap(sink.finish(bytes, statistics));
        } catch (const FailureException& error) {
            state.fail(error.failure);
        } catch (const std::exception& error) {
            state.fail(runtimeFailure(error.what()));
        }
    }
    if (state.error) {
        source.cancel();
        sink.cancel();
        return *state.error;
    }
    return report;
}

} // namespace nustar::star2711b
