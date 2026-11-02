#include "mapping_interface.hpp"

#include <utility>

namespace nustar::star2711b {
namespace {

Failure mappingError(ErrorKind kind, std::string message,
                     std::optional<FragmentId> fragment = {}) {
    return Failure{kind, Stage::Plan, fragment, 1, std::move(message)};
}

class PairedMapperWorkspace final : public StageWorkspace {
public:
    PairedMapperWorkspace(std::string semantic_id,
                          std::unique_ptr<StageWorkspace> seed,
                          std::unique_ptr<StageWorkspace> window,
                          std::unique_ptr<StageWorkspace> stitch)
        : identity_{"star-2.7.11b-paired-candidate-v1", 1,
                    std::move(semantic_id)},
          seed(std::move(seed)), window(std::move(window)),
          stitch(std::move(stitch)) {}
    const WorkspaceIdentity& identity() const noexcept override { return identity_; }
    void reset() noexcept override {
        seed->reset();
        window->reset();
        stitch->reset();
    }
    std::unique_ptr<StageWorkspace> seed, window, stitch;
    SeedBatch seeds;
    WindowBatch windows;
private:
    WorkspaceIdentity identity_;
};

class PairedCandidateMapper final : public IPairedCandidateMapper {
public:
    PairedCandidateMapper(ReferenceLease reference,
        std::shared_ptr<const ISeedFinder> seed,
        std::shared_ptr<const IWindowBuilder> window,
        std::shared_ptr<const IPairedStitcher> stitch,
        PairedCandidateMappingConfig config)
        : reference_(std::move(reference)), seed_(std::move(seed)),
          window_(std::move(window)), stitch_(std::move(stitch)),
          config_(std::move(config)),
          semantic_id_(reference_->view().semantic_id) {}

    std::string_view backendId() const noexcept override {
        return "star-2.7.11b-paired-candidate-v1";
    }

    Result<std::unique_ptr<StageWorkspace>> createWorkspace(
        const WorkspaceBudget& budget) const override {
        const auto view = reference_->view();
        auto seed = seed_->createWorkspace(config_.limits, budget);
        if (std::holds_alternative<Failure>(seed))
            return std::get<Failure>(std::move(seed));
        auto window = window_->createWorkspace(
            config_.window, view, config_.limits, budget);
        if (std::holds_alternative<Failure>(window))
            return std::get<Failure>(std::move(window));
        auto stitch = stitch_->createWorkspace(
            config_.stitch, config_.score, view, config_.limits, budget);
        if (std::holds_alternative<Failure>(stitch))
            return std::get<Failure>(std::move(stitch));
        return std::unique_ptr<StageWorkspace>(std::make_unique<PairedMapperWorkspace>(
            semantic_id_, std::get<std::unique_ptr<StageWorkspace>>(std::move(seed)),
            std::get<std::unique_ptr<StageWorkspace>>(std::move(window)),
            std::get<std::unique_ptr<StageWorkspace>>(std::move(stitch))));
    }

    Status mapPair(const PairedMappingRead& typed, StageWorkspace& workspace,
                   CandidateSet& out) const override {
        auto* stages = static_cast<PairedMapperWorkspace*>(&workspace);
        const auto& read = typed.view();
        stages->reset();
        auto& seeds = stages->seeds;
        if (auto status = seed_->find(read, *stages->seed, seeds);
            std::holds_alternative<Failure>(status)) return status;
        if (seeds.markers.normal_early_stop) {
            out = CandidateSet{};
            out.generation = reference_->view().generation;
            out.markers = seeds.markers;
            return std::monostate{};
        }
        auto& windows = stages->windows;
        const auto view = reference_->view();
        if (auto status = window_->build(read, seeds, view, config_.window,
                                          *stages->window, windows);
            std::holds_alternative<Failure>(status)) return status;
        if (windows.markers.normal_early_stop) {
            out = CandidateSet{};
            out.generation = view.generation;
            out.markers = windows.markers;
            return std::monostate{};
        }
        return stitch_->stitchPair(typed, windows, view, config_.stitch,
                                   config_.score, *stages->stitch, out);
    }
private:
    ReferenceLease reference_;
    std::shared_ptr<const ISeedFinder> seed_;
    std::shared_ptr<const IWindowBuilder> window_;
    std::shared_ptr<const IPairedStitcher> stitch_;
    PairedCandidateMappingConfig config_;
    std::string semantic_id_;
};

} // namespace

Result<std::shared_ptr<const IPairedCandidateMapper>> composePairedCandidateMapper(
    ReferenceLease reference, std::shared_ptr<const ISeedFinder> seed,
    std::shared_ptr<const IWindowBuilder> window,
    std::shared_ptr<const IPairedStitcher> stitch,
    const PairedCandidateMappingConfig& config) {
    if (!reference || !seed || !window || !stitch)
        return mappingError(ErrorKind::Config, "paired candidate mapper needs all kernels");
    if (auto status = window->validateSeedContract(seed->outputContract(), config.window);
        std::holds_alternative<Failure>(status))
        return std::get<Failure>(std::move(status));
    return std::shared_ptr<const IPairedCandidateMapper>(
        std::make_shared<PairedCandidateMapper>(std::move(reference),
            std::move(seed), std::move(window), std::move(stitch), config));
}

} // namespace nustar::star2711b
