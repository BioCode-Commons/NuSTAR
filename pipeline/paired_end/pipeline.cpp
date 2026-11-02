#include "pipeline.hpp"
#include "../../statistics/interface.hpp"

#include <array>
#include <utility>

namespace nustar::star2711b {
namespace {

Failure pipelineError(std::string message,
    std::optional<FragmentId> fragment = {}) {
    return Failure{ErrorKind::Config, Stage::Plan, fragment, 1,
                   std::move(message)};
}

class PairedPipelineWorkspace final : public StageWorkspace {
public:
    PairedPipelineWorkspace(std::string semantic_id,
        std::unique_ptr<StageWorkspace> mapper_workspace,
        std::unique_ptr<StageWorkspace> sam_workspace,
        std::unique_ptr<StageWorkspace> overlap_workspace)
        : identity_{"star-2.7.11b-paired-pipeline-v1", 1,
                    std::move(semantic_id)},
          mapper(std::move(mapper_workspace)), sam(std::move(sam_workspace)),
          overlap(std::move(overlap_workspace)) {}
    const WorkspaceIdentity& identity() const noexcept override {
        return identity_;
    }
    void reset() noexcept override {
        mapper->reset();
        sam->reset();
        if (overlap) overlap->reset();
    }
    std::unique_ptr<StageWorkspace> mapper, sam, overlap;
private:
    WorkspaceIdentity identity_;
};

class PairedReadPipeline final : public IPairedReadPipeline {
public:
    PairedReadPipeline(ReferenceLease reference, PairedKernelBundle bundle,
        PairedPipelineConfig config)
        : reference_(std::move(reference)), bundle_(std::move(bundle)),
          config_(std::move(config)),
          semantic_id_(reference_->view().semantic_id) {}

    std::string_view backendId() const noexcept override {
        return "star-2.7.11b-paired-pipeline-v1";
    }

    Result<std::unique_ptr<StageWorkspace>> createWorkspace(
        const WorkspaceBudget& budget) const override {
        auto mapper = bundle_.mapper->createWorkspace(budget);
        if (const auto* failure = std::get_if<Failure>(&mapper)) return *failure;
        auto sam = bundle_.sam->createWorkspace(config_.sam,
                                                config_.limits, budget);
        if (const auto* failure = std::get_if<Failure>(&sam)) return *failure;
        std::unique_ptr<StageWorkspace> overlap;
        if (config_.overlap.peOverlapNbasesMin != 0) {
            auto made = bundle_.overlap->createWorkspace(config_.overlap,
                *bundle_.merged_single_mapper, config_.limits, budget);
            if (const auto* failure = std::get_if<Failure>(&made)) return *failure;
            overlap = std::get<std::unique_ptr<StageWorkspace>>(std::move(made));
        }
        return std::unique_ptr<StageWorkspace>(
            std::make_unique<PairedPipelineWorkspace>(semantic_id_,
                std::get<std::unique_ptr<StageWorkspace>>(std::move(mapper)),
                std::get<std::unique_ptr<StageWorkspace>>(std::move(sam)),
                std::move(overlap)));
    }

    Result<FragmentEffects> computePair(PreparedReadPair&& read,
        StageWorkspace& workspace, LegacyRandomState& random) const override {
        auto computation = tracePair(std::move(read), workspace, random);
        if (const auto* failure = std::get_if<Failure>(&computation))
            return *failure;
        return std::move(std::get<PairedReadComputation>(computation).effects);
    }

    Result<PairedReadComputation> tracePair(PreparedReadPair&& read,
        StageWorkspace& workspace, LegacyRandomState& random) const override {
        auto* stages = static_cast<PairedPipelineWorkspace*>(&workspace);
        stages->reset();
        const auto& mapping = read.mapping.view();
        CandidateSet candidates{};
        if (auto status = bundle_.mapper->mapPair(read.mapping,
                *stages->mapper, candidates);
            std::holds_alternative<Failure>(status))
            return std::get<Failure>(std::move(status));
        OverlapState overlap{};
        if (config_.overlap.peOverlapNbasesMin != 0) {
            if (auto status = bundle_.overlap->refinePair(read.mapping,
                    reference_->view(), config_.overlap, config_.score,
                    *bundle_.merged_single_mapper, *stages->overlap,
                    candidates, overlap);
                std::holds_alternative<Failure>(status))
                return std::get<Failure>(std::move(status));
        }
        SelectionResult selection{};
        if (auto status = bundle_.selection->select(mapping, candidates,
                config_.select, random, selection);
            std::holds_alternative<Failure>(status))
            return std::get<Failure>(std::move(status));
        auto filtered = filterRead(mapping, candidates, selection, config_.filter);
        if (const auto* failure = std::get_if<Failure>(&filtered)) return *failure;
        const auto decision = std::get<FilterDecision>(std::move(filtered));
        auto junctions = collectJunctions(mapping, candidates, selection,
                                           decision, config_.sj);
        if (const auto* failure = std::get_if<Failure>(&junctions)) return *failure;
        auto output = bundle_.sam->projectPair(read, candidates, selection,
            decision, reference_->view(), config_.sam, config_.select,
            config_.stitch, *stages->sam);
        if (const auto* failure = std::get_if<Failure>(&output)) return *failure;
        auto statistics = deriveStatistics(mapping, candidates, selection,
            decision, std::get<ProjectedOutput>(output).disposition);
        if (const auto* failure = std::get_if<Failure>(&statistics)) return *failure;
        const std::array<std::string_view, 2> qualities{
            std::string_view(read.quality_for_output[0]).substr(mapping.mates[0].clipping.total_5p,mapping.mates[0].clipped_length),
            std::string_view(read.quality_for_output[1]).substr(mapping.mates[1].clipping.total_5p,mapping.mates[1].clipped_length)};
        FragmentEffects effects{read.original.id,
            std::get<FragmentJunctions>(std::move(junctions)),
            std::get<ProjectedOutput>(std::move(output)),
            std::get<StatsDelta>(std::move(statistics)),
            deriveQualityHistogram(qualities)};
        return PairedReadComputation{std::move(read), std::move(candidates),
            overlap, std::move(selection), decision, std::move(effects)};
    }

private:
    ReferenceLease reference_;
    PairedKernelBundle bundle_;
    PairedPipelineConfig config_;
    std::string semantic_id_;
};

} // namespace

Result<std::shared_ptr<const IPairedReadPipeline>> composePairedPipeline(
    ReferenceLease reference, PairedKernelBundle bundle,
    const PairedPipelineConfig& config) {
    if (!reference || !bundle.mapper || !bundle.selection || !bundle.sam)
        return pipelineError("paired pipeline needs reference, mapper, selection and SAM kernels");
    if (!reference->view().metadata || !reference->view().sequence)
        return pipelineError("paired pipeline needs a complete reference view");
    if (config.overlap.peOverlapNbasesMin != 0 &&
        (!bundle.overlap || !bundle.merged_single_mapper))
        return pipelineError("paired overlap needs refiner and merged-read mapper");
    return std::shared_ptr<const IPairedReadPipeline>(
        std::make_shared<PairedReadPipeline>(std::move(reference),
                                             std::move(bundle), config));
}

} // namespace nustar::star2711b
