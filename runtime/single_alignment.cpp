#include "single_alignment.hpp"
#include "../output/sink.hpp"
#include "../pipeline/single_end/mapping_interface.hpp"
#include "../pipeline/single_end/pipeline.hpp"
#include "../reads/barcode_source.hpp"
#include "../reads/fastq_prepare.hpp"
#include "../reads/limited_source.hpp"
#include "../sam/star_single.hpp"
#include "../seed_search/star_search.hpp"
#include "../selection/star_select.hpp"
#include "../solo/session.hpp"
#include "../splice_junction/star_junction.hpp"
#include "../splice_junction/table_writer.hpp"
#include "../statistics/interface.hpp"
#include "../stitching/single_stitch.hpp"
#include "../window/star_window.hpp"
#include "bam_output.hpp"
#include "defaults.hpp"
#include "input_tags.hpp"
#include "alignment_scheduler.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>

namespace nustar::star2711b {
namespace {

// Workers update these counters on every read. Do not let adjacent workers
// invalidate each other's cache line merely because profiling is enabled.
struct alignas(64) WorkerTiming {
    SingleSamReport values;
};

template <class T> T checked(Result<T> value) {
    if (std::holds_alternative<Failure>(value))
        throw std::runtime_error(std::get<Failure>(value).message);
    return std::get<T>(std::move(value));
}

void checked(Status value) {
    if (std::holds_alternative<Failure>(value))
        throw std::runtime_error(std::get<Failure>(value).message);
}

} // namespace

Result<SingleSamReport> runSingleSam(const SingleSamOptions& options) {
    try {
        if (options.threads == 0 || options.threads > 256)
            throw std::runtime_error("invalid worker count");
        const auto job_started = std::chrono::system_clock::now();
        const bool profile = options.profile_stages;
        using Clock = std::chrono::steady_clock;
        const auto elapsed = [](Clock::time_point start) {
            return std::chrono::duration<double>(Clock::now() - start).count();
        };
        SingleSamReport report{};
        auto initialization_started = Clock::now();
        auto reference = options.run.resources.reference
                             ? options.run.resources.reference
                             : checked(loadGenomeReference(options.genome_dir));
        report.reference_load_seconds = elapsed(initialization_started);
        SuffixArrayIndexBackend backend;
        SearchIndexRequest request{"nustar-sa", options.genome_dir / "index.sa", {}};
        initialization_started = Clock::now();
        auto index = options.run.resources.index ? options.run.resources.index
                                                 : checked(backend.load(request, reference));
        report.index_load_seconds = elapsed(initialization_started);
        const auto metadata = *reference->view().metadata;
        std::shared_ptr<solo::Session> solo_session;
        if (options.run.solo_config)
            solo_session = std::make_shared<solo::Session>(
                *options.run.solo_config, options.genome_dir,
                options.outFileNamePrefix.value_or(options.output_sam.string() + "."),
                options.threads, 1, metadata);
        const SeedDeliveryContract delivery{
            SeedContract::Star2711bEquivalent,   HitOrder::Star2711b, true, true, true, 1,
            std::string(star_seed_search_domain)};
        const auto seed_config = defaultStarSeedConfig();
        auto seed_plan =
            std::make_shared<SeedSearchPlan>(SeedSearchPlan{delivery, {}, seed_config, {}});
        StarSeedAlgorithmFactory seed_factory;
        initialization_started = Clock::now();
        auto finder = checked(seed_factory.bind(index, reference, seed_plan));
        report.seed_bind_seconds = elapsed(initialization_started);
        initialization_started = Clock::now();
        auto builder = checked(makeStarWindowBuilder(reference));
        auto stitcher = checked(makeStarSingleStitcher(reference));
        const auto& storage = reference->storage();
        const auto win_bits = std::min<std::uint32_t>(16, storage.chromosome_bin_bits);
        const WindowConfig window_config{
            50,
            win_bits,
            storage.chromosome_bin_bits - win_bits,
            (storage.forward.size() + reference->synthetic().size()) / (1ULL << win_bits) + 1,
            9,
            4,
            50,
            10,
            10000,
            0,
            0};
        const auto stitch_config = defaultStarSingleStitchConfig();
        const auto score = defaultStarScoreConfig();
        const CompatibilityLimits limits{650, 1024, 20, 10000, 255};
        const WorkspaceBudget budget{std::numeric_limits<std::uint64_t>::max(), 128};
        const SingleCandidateMappingConfig candidate_config{window_config, stitch_config, score,
                                                            limits};
        auto mapper = checked(
            composeSingleCandidateMapper(reference, finder, builder, stitcher, candidate_config));
        auto selector = makeStarAlignmentSelector();
        const auto selection = defaultStarSelectionConfig();
        const auto filter = defaultStarReadFilterConfig();
        auto sam = defaultStarSamConfig();
        if (options.run.solo_config) {
            sam.read_group_by_file = options.run.solo_config->read_groups;
            sam.read_group_lines = options.run.solo_config->read_group_lines;
        }
        sam.enabled = options.run.sam_enabled;
        if (!options.run.output_attributes.empty())
            sam.outSAMattributes = options.run.output_attributes;
        sam.unmapped_within = options.outSAMunmappedWithin;
        auto projector = makeStarSingleSamProjector();
        const auto sj_config = defaultStarSjOutputConfig();
        const SinglePipelineConfig pipeline_config{selection, filter,        sj_config,
                                                   sam,       stitch_config, limits};
        // Composition validates selection once, also for the direct fast loop below.
        auto pipeline = checked(composeSinglePipeline(
            reference, SingleKernelBundle{mapper, selector, projector}, pipeline_config));
        // Like STAR's genomeMain.freeMemory(), drop mapping-only leases before
        // Solo reduction and BAM sorting. Metadata remains available to writers.
        auto release_mapping = [&] {
            pipeline.reset();
            mapper.reset();
            finder.reset();
            builder.reset();
            stitcher.reset();
            index.reset();
            reference.reset();
        };
        sam.command_line_full = options.command_line;
        sam.command_line_user = sam.command_line_full;
        const BaselineIdentity baseline{"2.7.11b", {}, {}, {}, {}, {}};
        const auto sam_header = checked(encodeSamHeader(metadata, sam, baseline));
        report.pipeline_setup_seconds = elapsed(initialization_started);
        if (!profile) {
            InputSpec input_spec{Layout::SE, {}, {}, {}};
            input_spec.groups = options.run.input_groups;
            if (input_spec.groups.empty())
                input_spec.groups.push_back(InputGroup{{options.reads_fastq, {}}});
            input_spec.capture_aux = options.run.bam_input;
            const auto& input_backend = options.run.bam_input ? bamSourceBackend()
                                                        : fastqSourceBackend();
            if (options.io_backend != "pread" && options.io_backend != "mmap" && options.io_backend != "buffered")
                throw std::runtime_error("unknown FASTQ IO backend");
            auto rna_source =
                limitReads(checked(input_backend.open(input_spec)), options.run.read_limit);
            auto source =
                checked(withBarcode(std::move(rna_source), input_backend, options.run.barcode_fastq,
                                    options.run.read_limit == UINT64_MAX));
            const auto mapping_started = std::chrono::system_clock::now();
            auto sink = checked(
                openAlignmentSink(OutputFiles{options.output_sam, job_started, mapping_started,
                                              options.outFileNamePrefix, options.run.sam_enabled},
                                  options.run.bam, metadata, solo_session));
            const auto coarse_started = Clock::now();
            auto aligned = checked(runParallelAlignment(
                *source, *sink, sam_header, sj_config, metadata, options.threads, 1ULL << 20,
                [pipeline, budget, solo_session]() -> Result<std::unique_ptr<StageWorkspace>> {
                    auto workspace = pipeline->createWorkspace(budget);
                    if (!solo_session || std::holds_alternative<Failure>(workspace))
                        return workspace;
                    return std::unique_ptr<StageWorkspace>(std::make_unique<solo::Workspace>(
                        std::get<std::unique_ptr<StageWorkspace>>(std::move(workspace)),
                        solo_session->createWorker()));
                },
                [pipeline, solo_session, clips = options.run.clipping,
                 solo_config = options.run.solo_config, keep = options.run.input_attributes,
                 bam_output =
                     options.run.bam.enabled()](RawFragment&& raw, StageWorkspace& workspace,
                                                LegacyRandomState& rng) -> Result<FragmentEffects> {
                    // The bound source supplies complete SE records with qualities.
                    auto& mate = *raw.mates[0];
                    auto read = prepareSingleFastq(std::move(mate.header), std::move(mate.sequence),
                                                   std::move(*mate.quality), raw.id.ordinal,
                                                   clips[0], raw.id.file_group);
                    read.original.legacy_read_filter = raw.legacy_read_filter;
                    if (!solo_session) {
                        auto result = pipeline->computeSingle(std::move(read), workspace, rng);
                        if (auto* effects = std::get_if<FragmentEffects>(&result))
                            appendInputTags(effects->output, raw, keep, bam_output);
                        return result;
                    }
                    auto& stages = static_cast<solo::Workspace&>(workspace);
                    auto result = pipeline->traceSingle(std::move(read), *stages.mapping, rng);
                    if (auto failure = std::get_if<Failure>(&result))
                        return *failure;
                    auto& computation = std::get<SingleReadComputation>(result);
                    stages.counting->record(
                        raw.barcode                 ? &*raw.barcode
                        : solo_config->barcode_mate ? &computation.read.original.mate
                                                    : nullptr,
                        computation.read.mapping.view(), computation.candidates,
                        computation.selection, computation.filter, computation.effects.statistics,
                        computation.effects.quality, &computation.effects.output,
                        mate.header_extra);
                    appendInputTags(computation.effects.output, raw, keep, bam_output);
                    return std::move(computation.effects);
                }));
            source.reset();
            release_mapping();
            if (solo_session)
                solo_session->finish();
            checked(sink->finalize());
            report.reads = aligned.fragments;
            report.records = aligned.sam_records;
            report.mapping_wall_seconds = elapsed(coarse_started);
            return report;
        }
        InputSpec input_spec{Layout::SE, {}, {}, {}};
        input_spec.groups = options.run.input_groups;
        if (input_spec.groups.empty())
            input_spec.groups.push_back(InputGroup{{options.reads_fastq, {}}});
        if (options.io_backend != "pread" && options.io_backend != "mmap" && options.io_backend != "buffered")
            throw std::runtime_error("unknown FASTQ IO backend");
        input_spec.capture_aux = options.run.bam_input;
        const auto& input_backend = options.run.bam_input ? bamSourceBackend()
                                                        : fastqSourceBackend();
        auto rna_source =
            limitReads(checked(input_backend.open(input_spec)), options.run.read_limit);
        auto source =
            checked(withBarcode(std::move(rna_source), input_backend, options.run.barcode_fastq,
                                options.run.read_limit == UINT64_MAX));
        const auto mapping_started = std::chrono::system_clock::now();
        auto sink = checked(
            openAlignmentSink(OutputFiles{options.output_sam, job_started, mapping_started,
                                          options.outFileNamePrefix, options.run.sam_enabled},
                              options.run.bam, metadata, solo_session));
        checked(sink->begin(sam_header));
        const auto mapping_wall_started = Clock::now();
        const bool concurrent_output = sink->concurrentWrites();
        checked(source->prepareParallel(1ULL << 20, options.threads));
        const bool concurrent_input = source->concurrentReads();
        std::vector<WorkerAlignmentReduction> reductions(options.threads);
        // Each worker owns its workspaces, RNG and reductions. Completed
        // batches submit immediately; only non-concurrent sinks need a lock.
        std::vector<WorkerTiming> timings(options.threads);
        std::mutex input_mutex, output_mutex, failure_mutex;
        std::atomic<bool> stopped{false};
        std::string error;
        const auto fail = [&](const std::exception& failure) {
            {
                std::lock_guard lock(failure_mutex);
                if (error.empty())
                    error = failure.what();
                stopped.store(true);
            }
        };
        std::vector<std::jthread> workers;
        try {
            workers.reserve(options.threads);
            for (std::uint32_t worker = 0; worker < options.threads; ++worker) {
                workers.emplace_back([&, worker] {
                    try {
                        auto counting = solo_session ? solo_session->createWorker() : nullptr;
                        auto seed_workspace = checked(finder->createWorkspace(limits, budget));
                        auto window_workspace = checked(builder->createWorkspace(
                            window_config, reference->view(), limits, budget));
                        auto stitch_workspace = checked(stitcher->createWorkspace(
                            stitch_config, score, reference->view(), limits, budget));
                        auto sam_workspace =
                            checked(projector->createWorkspace(sam, limits, budget));
                        LegacyRandomState random{};
                        random.engine.seed(777);
                        auto& timing = timings[worker].values;
                        SeedBatch seeds{};
                        WindowBatch windows{};
                        while (!stopped.load()) {
                            std::optional<ReadBatch> task;
                            if (concurrent_input) {
                                const auto reading = Clock::now();
                                task = checked(source->next(1ULL << 20));
                                timing.input_read_seconds += elapsed(reading);
                            } else {
                                const auto waiting = Clock::now();
                                std::lock_guard lock(input_mutex);
                                timing.input_wait_seconds += elapsed(waiting);
                                if (stopped.load())
                                    break;
                                const auto reading = Clock::now();
                                task = checked(source->next(1ULL << 20));
                                timing.input_read_seconds += elapsed(reading);
                            }
                            if (!task)
                                break;
                            const auto previous_records = timing.records;
                            OutputBatch output{task->id, {}};
                            std::vector<FragmentJunctions> batch_junctions;
                            std::vector<StatsDelta> batch_statistics;
                            output.fragments.reserve(task->fragments.size());
                            batch_junctions.reserve(task->fragments.size());
                            batch_statistics.reserve(task->fragments.size());
                            for (auto& raw : task->fragments) {
                                if (stopped.load())
                                    return;
                                // Layout and qualities are established by the bound source.
                                auto& mate = *raw.mates[0];
                                const auto preparing = Clock::now();
                                auto read = prepareSingleFastq(
                                    std::move(mate.header), std::move(mate.sequence),
                                    std::move(*mate.quality), raw.id.ordinal,
                                    options.run.clipping[0], raw.id.file_group);
                                read.original.legacy_read_filter = raw.legacy_read_filter;
                                timing.prepare_seconds += elapsed(preparing);
                                CandidateSet candidates{};
                                auto started = Clock::now();
                                checked(finder->find(read.mapping.view(), *seed_workspace, seeds));
                                timing.seed_seconds += elapsed(started);
                                if (seeds.markers.normal_early_stop) {
                                    candidates.generation = reference->view().generation;
                                    candidates.markers = seeds.markers;
                                } else {
                                    started = Clock::now();
                                    checked(builder->build(read.mapping.view(), seeds,
                                                           reference->view(), window_config,
                                                           *window_workspace, windows));
                                    timing.window_seconds += elapsed(started);
                                    if (windows.markers.normal_early_stop) {
                                        candidates.generation = windows.generation;
                                        candidates.markers = windows.markers;
                                    } else {
                                        started = Clock::now();
                                        checked(stitcher->stitchSingle(
                                            read.mapping, windows, reference->view(), stitch_config,
                                            score, *stitch_workspace, candidates));
                                        timing.stitch_seconds += elapsed(started);
                                    }
                                }
                                SelectionResult selected{};
                                started = Clock::now();
                                checked(selector->select(read.mapping.view(), candidates, selection,
                                                         random, selected));
                                const auto decision = checked(
                                    filterRead(read.mapping.view(), candidates, selected, filter));
                                batch_junctions.push_back(
                                    checked(collectJunctions(read.mapping.view(), candidates,
                                                             selected, decision, sj_config)));
                                timing.selection_seconds += elapsed(started);
                                started = Clock::now();
                                auto projected = checked(projector->projectSingle(
                                    read, candidates, selected, decision, reference->view(), sam,
                                    selection, stitch_config, *sam_workspace));
                                batch_statistics.push_back(checked(
                                    deriveStatistics(read.mapping.view(), candidates, selected,
                                                     decision, projected.disposition)));
                                if (counting) {
                                    const auto& mate = read.mapping.view().mates[0];
                                    const std::array<std::string_view, 1> quality{
                                        std::string_view(read.quality_for_output)
                                            .substr(mate.clipping.total_5p, mate.clipped_length)};
                                    counting->record(raw.barcode ? &*raw.barcode
                                                     : options.run.solo_config->barcode_mate
                                                         ? &read.original.mate
                                                         : nullptr,
                                                     read.mapping.view(), candidates, selected,
                                                     decision, batch_statistics.back(),
                                                     deriveQualityHistogram(quality), &projected,
                                                     raw.mates[0]->header_extra);
                                }
                                appendInputTags(projected, raw, options.run.input_attributes,
                                                options.run.bam.enabled());
                                output.fragments.push_back(
                                    checked(sink->encode(read.original.id, projected)));
                                ++timing.reads;
                                timing.records += projected.records.size();
                                timing.sam_seconds += elapsed(started);
                            }
                            const auto writing = Clock::now();
                            if (concurrent_output) {
                                checked(sink->submit(std::move(output)));
                            } else {
                                const auto waiting = Clock::now();
                                std::unique_lock lock(output_mutex);
                                timing.commit_wait_seconds += elapsed(waiting);
                                checked(sink->submit(std::move(output)));
                            }
                            reductions[worker].add(task->id, batch_junctions, batch_statistics,
                                                   timing.records - previous_records);
                            timing.sam_seconds += elapsed(writing);
                        }
                    } catch (const std::exception& failure) {
                        fail(failure);
                    }
                });
            }
        } catch (const std::exception& failure) {
            fail(failure);
        }
        for (auto& worker : workers)
            worker.join();
        if (stopped.load()) {
            source->cancel();
            sink->cancel();
            throw std::runtime_error(error);
        }
        checked(finishConcurrentAlignment(reductions, *sink, sj_config, metadata));
        source.reset();
        release_mapping();
        if (solo_session)
            solo_session->finish();
        checked(sink->finalize());
        report.mapping_wall_seconds = elapsed(mapping_wall_started);
        for (const auto& worker_timing : timings) {
            const auto& timing = worker_timing.values;
            report.reads += timing.reads;
            report.records += timing.records;
            report.seed_seconds += timing.seed_seconds;
            report.window_seconds += timing.window_seconds;
            report.stitch_seconds += timing.stitch_seconds;
            report.selection_seconds += timing.selection_seconds;
            report.sam_seconds += timing.sam_seconds;
            report.input_wait_seconds += timing.input_wait_seconds;
            report.commit_wait_seconds += timing.commit_wait_seconds;
            report.prepare_seconds += timing.prepare_seconds;
            report.input_read_seconds += timing.input_read_seconds;
        }
        return report;
    } catch (const std::exception& error) {
        return Failure{ErrorKind::Internal, Stage::Finish, {}, 1, error.what()};
    }
}

} // namespace nustar::star2711b
