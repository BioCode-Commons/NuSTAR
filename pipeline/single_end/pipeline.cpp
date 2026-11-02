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

class SinglePipelineWorkspace final : public StageWorkspace {
public:
    SinglePipelineWorkspace(std::string semantic_id,
        std::unique_ptr<StageWorkspace> mapper_workspace,
        std::unique_ptr<StageWorkspace> sam_workspace)
        : identity_{"star-2.7.11b-single-pipeline-v1", 1,
                    std::move(semantic_id)},
          mapper(std::move(mapper_workspace)), sam(std::move(sam_workspace)) {}
    const WorkspaceIdentity& identity() const noexcept override {
        return identity_;
    }
    void reset() noexcept override {
        mapper->reset();
        sam->reset();
    }
    std::unique_ptr<StageWorkspace> mapper, sam;
private:
    WorkspaceIdentity identity_;
};

class SingleReadPipeline final : public ISingleReadPipeline {
public:
    SingleReadPipeline(ReferenceLease reference, SingleKernelBundle bundle,
        SinglePipelineConfig config)
        : reference_(std::move(reference)), bundle_(std::move(bundle)),
          config_(std::move(config)),
          semantic_id_(reference_->view().semantic_id) {}

    std::string_view backendId() const noexcept override {
        return "star-2.7.11b-single-pipeline-v1";
    }

    Result<std::unique_ptr<StageWorkspace>> createWorkspace(
        const WorkspaceBudget& budget) const override {
        auto mapper = bundle_.mapper->createWorkspace(budget);
        if (const auto* failure = std::get_if<Failure>(&mapper)) return *failure;
        auto sam = bundle_.sam->createWorkspace(config_.sam,
                                                config_.limits, budget);
        if (const auto* failure = std::get_if<Failure>(&sam)) return *failure;
        return std::unique_ptr<StageWorkspace>(
            std::make_unique<SinglePipelineWorkspace>(semantic_id_,
                std::get<std::unique_ptr<StageWorkspace>>(std::move(mapper)),
                std::get<std::unique_ptr<StageWorkspace>>(std::move(sam))));
    }

    Result<FragmentEffects> computeSingle(PreparedSingleRead&& read,
        StageWorkspace& workspace, LegacyRandomState& random) const override {
        auto computation = traceSingle(std::move(read), workspace, random);
        if (const auto* failure = std::get_if<Failure>(&computation))
            return *failure;
        return std::move(std::get<SingleReadComputation>(computation).effects);
    }

    Result<SingleReadComputation> traceSingle(PreparedSingleRead&& read,
        StageWorkspace& workspace, LegacyRandomState& random) const override {
        auto* stages = static_cast<SinglePipelineWorkspace*>(&workspace);
        stages->reset();
        const auto& mapping = read.mapping.view();
        CandidateSet candidates{};
        if (auto status = bundle_.mapper->mapSingle(read.mapping,
                *stages->mapper, candidates);
            std::holds_alternative<Failure>(status))
            return std::get<Failure>(std::move(status));
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
        auto output = bundle_.sam->projectSingle(read, candidates, selection,
            decision, reference_->view(), config_.sam, config_.select,
            config_.stitch, *stages->sam);
        if (const auto* failure = std::get_if<Failure>(&output)) return *failure;
        auto statistics = deriveStatistics(mapping, candidates, selection,
            decision, std::get<ProjectedOutput>(output).disposition);
        if (const auto* failure = std::get_if<Failure>(&statistics)) return *failure;
        const std::array<std::string_view, 1> qualities{std::string_view(read.quality_for_output).substr(
            mapping.mates[0].clipping.total_5p,mapping.mates[0].clipped_length)};
        FragmentEffects effects{read.original.id,
            std::get<FragmentJunctions>(std::move(junctions)),
            std::get<ProjectedOutput>(std::move(output)),
            std::get<StatsDelta>(std::move(statistics)),
            deriveQualityHistogram(qualities)};
        return SingleReadComputation{std::move(read), std::move(candidates),
            std::move(selection), decision, std::move(effects)};
    }

private:
    ReferenceLease reference_;
    SingleKernelBundle bundle_;
    SinglePipelineConfig config_;
    std::string semantic_id_;
};

} // namespace

Result<std::shared_ptr<const ISingleReadPipeline>> composeSinglePipeline(
    ReferenceLease reference, SingleKernelBundle bundle,
    const SinglePipelineConfig& config) {
    if (!reference || !bundle.mapper || !bundle.selection || !bundle.sam)
        return pipelineError("single pipeline needs reference, mapper, selection and SAM kernels");
    if (!reference->view().metadata || !reference->view().sequence)
        return pipelineError("single pipeline needs a complete reference view");
    return std::shared_ptr<const ISingleReadPipeline>(
        std::make_shared<SingleReadPipeline>(std::move(reference),
                                             std::move(bundle), config));
}

} // namespace nustar::star2711b
