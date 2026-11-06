#include "paired_alignment.hpp"
#include "../output/sink.hpp"
#include "../pipeline/paired_end/mapping_interface.hpp"
#include "../pipeline/paired_end/pipeline.hpp"
#include "../reads/barcode_source.hpp"
#include "../reads/fastq_prepare.hpp"
#include "../reads/limited_source.hpp"
#include "../sam/star_paired.hpp"
#include "../seed_search/star_search.hpp"
#include "../selection/star_select.hpp"
#include "../solo/session.hpp"
#include "../splice_junction/star_junction.hpp"
#include "../statistics/interface.hpp"
#include "../stitching/paired_stitch.hpp"
#include "../stitching/single_stitch.hpp"
#include "../window/star_window.hpp"
#include "bam_output.hpp"
#include "defaults.hpp"
#include "input_tags.hpp"
#include "ordered_alignment.hpp"

#include <algorithm>
#include <chrono>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace nustar::star2711b {
namespace {

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

Result<PairedSamReport> runPairedSam(const PairedSamOptions& options) {
    const auto initialization_started = std::chrono::steady_clock::now();
    try {
        const auto job_started = std::chrono::system_clock::now();
        auto reference = options.run.resources.reference
                             ? options.run.resources.reference
                             : checked(loadGenomeReference(options.genome_dir));
        SuffixArrayIndexBackend backend;
        SearchIndexRequest request{"nustar-sa", options.genome_dir / "index.sa", {}};
        auto index = options.run.resources.index ? options.run.resources.index
                                                 : checked(backend.load(request, reference));
        const auto metadata = *reference->view().metadata;
        std::shared_ptr<solo::Session> solo_session;
        if (options.run.solo_config)
            solo_session = std::make_shared<solo::Session>(
                *options.run.solo_config, options.genome_dir,
                options.outFileNamePrefix.value_or(options.output_sam.string() + "."),
                options.threads, 2, metadata);
        const SeedDeliveryContract delivery{
            SeedContract::Star2711bEquivalent,   HitOrder::Star2711b, true, true, true, 1,
            std::string(star_seed_search_domain)};
        const auto seed_config = defaultStarSeedConfig();
        auto seed_plan =
            std::make_shared<SeedSearchPlan>(SeedSearchPlan{delivery, {}, seed_config, {}});
        StarSeedAlgorithmFactory seed_factory;
        auto finder = checked(seed_factory.bind(index, reference, seed_plan));
        auto builder = checked(makeStarWindowBuilder(reference));
        auto stitcher = checked(makeStarPairedStitcher(reference));
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
        const auto stitch_config = defaultStarPairedStitchConfig();
        const auto score = defaultStarScoreConfig();
        const CompatibilityLimits limits{650, 1024, 20, 10000, 255};
        const WorkspaceBudget budget{std::numeric_limits<std::uint64_t>::max(), 128};
        const PairedCandidateMappingConfig candidate_config{window_config, stitch_config, score,
                                                            limits};
        auto mapper = checked(
            composePairedCandidateMapper(reference, finder, builder, stitcher, candidate_config));
        std::shared_ptr<const ISingleCandidateMapper> merged_mapper;
        std::shared_ptr<const IPairedOverlapRefiner> overlap;
        if (options.run.overlap.peOverlapNbasesMin) {
            merged_mapper=checked(composeSingleCandidateMapper(reference,finder,builder,
                checked(makeStarSingleStitcher(reference)),
                SingleCandidateMappingConfig{window_config,defaultStarSingleStitchConfig(),score,limits}));
            overlap=checked(makeStarPairedOverlapRefiner(reference));
        }
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
        sam.outSAMorder = options.outSAMorder;
        PairedSamConfig paired_sam{sam, false};
        auto projector = makeStarPairedSamProjector();
        const auto sj_config = defaultStarSjOutputConfig();
        const PairedPipelineConfig pipeline_config{selection,  filter,        sj_config,
                                                   paired_sam, stitch_config, options.run.overlap,
                                                   score,      limits};
        auto pipeline = checked(composePairedPipeline(
            reference, PairedKernelBundle{mapper, merged_mapper, overlap, selector, projector}, pipeline_config));
        // Like STAR's genomeMain.freeMemory(), drop mapping-only leases before
        // Solo reduction and BAM sorting. Metadata remains available to writers.
        auto release_mapping = [&] {
            pipeline.reset();
            mapper.reset();
            merged_mapper.reset();
            overlap.reset();
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
        InputSpec input_spec{Layout::PE, {}, {}, {}};
        input_spec.groups = options.run.input_groups;
        if (input_spec.groups.empty())
            input_spec.groups.push_back(InputGroup{{options.mate1_fastq, options.mate2_fastq}});
        input_spec.capture_aux = options.run.bam_input;
        const auto& input_backend = options.run.bam_input          ? bamSourceBackend()
                                    : options.io_backend == "mmap" ? memoryMappedFastqBackend()
                                                                   : bufferedFastqBackend();
        if (options.io_backend != "mmap" && options.io_backend != "buffered")
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
        const auto coarse_started = std::chrono::steady_clock::now();
        auto report = checked(runOrderedAlignment(
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
                if (raw.layout != Layout::PE || !raw.mates[0] || !raw.mates[1])
                    return Failure{ErrorKind::Input, Stage::Input, raw.id, 1,
                                   "invalid paired FASTQ fragment"};
                auto& first = *raw.mates[0];
                auto& second = *raw.mates[1];
                auto read = preparePairedFastq(
                    std::move(first.header), std::move(first.sequence), std::move(*first.quality),
                    std::move(second.header), std::move(second.sequence),
                    std::move(*second.quality), raw.id.ordinal, clips, raw.id.file_group);
                read.original.legacy_read_filter = raw.legacy_read_filter;
                if (!solo_session) {
                    auto result = pipeline->computePair(std::move(read), workspace, rng);
                    if (auto* effects = std::get_if<FragmentEffects>(&result))
                        appendInputTags(effects->output, raw, keep, bam_output);
                    return result;
                }
                auto& stages = static_cast<solo::Workspace&>(workspace);
                auto result = pipeline->tracePair(std::move(read), *stages.mapping, rng);
                if (auto failure = std::get_if<Failure>(&result))
                    return *failure;
                auto& computation = std::get<PairedReadComputation>(result);
                stages.counting->record(
                    raw.barcode ? &*raw.barcode
                    : solo_config->barcode_mate
                        ? &computation.read.original.mates[solo_config->barcode_mate - 1]
                        : nullptr,
                    computation.read.mapping.view(), computation.candidates, computation.selection,
                    computation.filter, computation.effects.statistics, computation.effects.quality,
                    &computation.effects.output, first.header_extra);
                appendInputTags(computation.effects.output, raw, keep, bam_output);
                return std::move(computation.effects);
            }));
        source.reset();
        release_mapping();
        if (solo_session)
            solo_session->finish();
        checked(sink->finalize());
        return PairedSamReport{
            report.fragments, report.sam_records,
            std::chrono::duration<double>(coarse_started - initialization_started).count(),
            std::chrono::duration<double>(std::chrono::steady_clock::now() - coarse_started)
                .count()};
    } catch (const std::exception& error) {
        return Failure{ErrorKind::Internal, Stage::Finish, {}, 1, error.what()};
    }
}

} // namespace nustar::star2711b
