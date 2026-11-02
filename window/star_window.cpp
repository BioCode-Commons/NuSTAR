#include "star_window.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <utility>

namespace nustar::star2711b {
namespace {

constexpr std::uint16_t no_window = std::numeric_limits<std::uint16_t>::max();

Failure windowError(ErrorKind kind, std::string message,
                    std::optional<FragmentId> fragment = {}) {
    return Failure{kind, Stage::Window, fragment, 1, std::move(message)};
}

class StarWindowWorkspace final : public StageWorkspace {
public:
    StarWindowWorkspace(std::string semantic_id, const WindowConfig& config,
                        const ReferenceMetadata& metadata)
        : identity_{"star-2.7.11b-window-v1", 1, std::move(semantic_id)} {
        for (auto& strand : bin_to_window)
            strand.resize(config.winBinN, no_window);
        chromosome_for_bin.reserve(config.winBinN);
        const auto& chromosomes = metadata.chromosomes;
        for (std::uint64_t bin = 0; bin < config.winBinN; ++bin) {
            const auto position = (bin >> config.winBinChrNbits)
                << (config.winBinChrNbits + config.winBinNbits);
            auto it = std::upper_bound(chromosomes.begin(), chromosomes.end(),
                position, [](std::uint64_t value, const Chromosome& chr) {
                    return value < chr.start.value;
                });
            chromosome_for_bin.push_back(it == chromosomes.begin() ?
                chromosomes.front().id : std::prev(it)->id);
        }
    }
    const WorkspaceIdentity& identity() const noexcept override { return identity_; }
    void reset() noexcept override {
        for (unsigned strand = 0; strand < 2; ++strand) {
            for (const auto bin : touched_bins[strand])
                bin_to_window[strand][bin] = no_window;
            touched_bins[strand].clear();
        }
    }
    void setBin(Orientation orientation, std::uint64_t bin, std::uint16_t window) {
        const auto strand = static_cast<unsigned>(orientation);
        auto& previous = bin_to_window[strand][bin];
        // Record before mutation so allocation failure leaves reset() safe.
        // Rewrites during merges/flank extension need no additional entry.
        if (previous == no_window && window != no_window)
            touched_bins[strand].push_back(bin);
        previous = window;
    }
    std::array<std::vector<std::uint16_t>, 2> bin_to_window;
    std::vector<ChromosomeId> chromosome_for_bin;
private:
    std::array<std::vector<std::uint64_t>, 2> touched_bins;
    WorkspaceIdentity identity_;
};

class StarWindowBuilder final : public IWindowBuilder {
public:
    explicit StarWindowBuilder(ReferenceLease reference)
        : reference_(std::move(reference)), semantic_id_(reference_->view().semantic_id) {}

    std::string_view backendId() const noexcept override {
        return "star-2.7.11b-window-v1";
    }

    Status validateSeedContract(const SeedDeliveryContract& seed,
                                const WindowConfig&) const override {
        if (seed.strategy != SeedContract::Star2711bEquivalent ||
            seed.hit_order != HitOrder::Star2711b ||
            !seed.exact_multiplicity || !seed.complete_occurrences ||
            !seed.known_junction_projection || seed.hit_schema_version != 1)
            return windowError(ErrorKind::Config, "window builder requires STAR ordered seed hits");
        return std::monostate{};
    }

    Result<std::unique_ptr<StageWorkspace>> createWorkspace(
        const WindowConfig& config, const ReferenceView& reference,
        const CompatibilityLimits&, const WorkspaceBudget&) const override {
        if (reference.semantic_id != semantic_id_ || !reference.metadata ||
            reference.metadata->chromosomes.empty() || config.winBinN == 0 ||
            config.winBinNbits >= 63 || config.winBinChrNbits >= 63 ||
            config.winBinNbits+config.winBinChrNbits >= 63 ||
            config.alignWindowsPerReadNmax == 0 ||
            config.alignWindowsPerReadNmax > no_window ||
            config.seedPerWindowNmax == 0)
            return windowError(ErrorKind::Config, "invalid STAR window reference or capacity");
        return std::unique_ptr<StageWorkspace>(
            std::make_unique<StarWindowWorkspace>(semantic_id_, config,
                *reference.metadata));
    }

    Status build(const MappingRead& read, const SeedBatch& seeds,
                 const ReferenceView& reference, const WindowConfig& config,
                 StageWorkspace& workspace, WindowBatch& out) const override {
        auto* scratch = static_cast<StarWindowWorkspace*>(&workspace);
        scratch->reset();
        out.generation = reference.generation;
        out.markers = seeds.markers;
        if (seeds.markers.normal_early_stop) {
            out.windows.clear();
            return std::monostate{};
        }
        auto& windows = out.windows;
        std::size_t window_count = 0;
        const auto& chromosomes = reference.metadata->chromosomes;


        // STAR's chrBin is indexed by genomeChrBinNbases. The chromosome
        // boundaries are bin-aligned, so its value follows the latest start.
        const auto chromosomeForBin = [&](std::uint64_t bin) -> ChromosomeId {
            return scratch->chromosome_for_bin[bin];
        };

        const auto createWindow = [&](std::uint64_t genome_start,
                                      Orientation orientation) -> bool {
            const auto bin = genome_start >> config.winBinNbits;
            auto& map = scratch->bin_to_window[static_cast<unsigned>(orientation)];
            if (map[bin] != no_window) return true;
            const auto chromosome = chromosomeForBin(bin);
            bool left = false, right = false;
            std::uint64_t left_bin = bin, right_bin = bin;
            std::uint16_t window = no_window, right_window = no_window;
            if (bin > 0) {
                const auto minimum = bin > config.winAnchorDistNbins ?
                    bin-config.winAnchorDistNbins : 0;
                for (std::uint64_t probe = bin; probe-- > minimum;) {
                    if (map[probe] != no_window) {
                        left = chromosomeForBin(probe).value == chromosome.value;
                        if (left) {
                            window = map[probe];
                            left_bin = windows[window].first_bin.value;
                            for (auto fill = probe+1; fill <= bin; ++fill)
                                scratch->setBin(orientation, fill, window);
                        }
                        break;
                    }
                }
            }
            if (bin+1 < config.winBinN) {
                const auto maximum = std::min(config.winBinN,
                    bin+config.winAnchorDistNbins+1);
                for (auto probe = bin+1; probe < maximum; ++probe) {
                    if (map[probe] != no_window) {
                        right = chromosomeForBin(probe).value == chromosome.value;
                        if (right) {
                            right_window = map[probe];
                            right_bin = windows[right_window].last_bin_inclusive.value;
                            if (!left) window = right_window;
                            for (auto fill = bin; fill <= right_bin; ++fill)
                                scratch->setBin(orientation, fill, window);
                        }
                        break;
                    }
                }
            }
            if (!left && !right) {
                const auto id = window_count++;
                if (id == windows.size()) windows.emplace_back();
                auto& slot = windows[id];
                slot.seeds.clear();
                slot = AlignmentWindow{WindowId{id}, chromosome,
                    orientation, BinId{bin}, BinId{bin}, true, std::move(slot.seeds), 0, {}};
                scratch->setBin(orientation, bin, static_cast<std::uint16_t>(id));
                // STAR drops the just-created window at the configured limit.
                if (window_count >= config.alignWindowsPerReadNmax) {
                    --window_count;
                    scratch->setBin(orientation, bin, no_window);
                    return false;
                }
            } else {
                auto& existing = windows[window];
                existing.first_bin = BinId{left_bin};
                existing.last_bin_inclusive = BinId{right_bin};
                if (left && right) windows[right_window].alive = false;
            }
            return true;
        };

        // Anchor pass: PC order, then SA rank order. A crossing SJ contributes
        // donor and acceptor anchors; a non-crossing synthetic occurrence does not.
        for (const auto& seed : seeds.seeds) {
            if (seed.multiplicity.value > config.winAnchorMultimapNmax) continue;
            auto status = seeds.occurrences->visitAnchors(seed.occurrences,
                [&](GenomeOffset start, Orientation orientation) {
                    return createWindow(start.value, orientation);
                });
            if (std::holds_alternative<Failure>(status)) return status;
        }

        windows.resize(window_count);
        for (auto& window : windows) {
            if (!window.alive) continue;
            auto bin = window.first_bin.value;
            for (std::uint64_t i = 0; i < config.winFlankNbins && bin > 0 &&
                 chromosomeForBin(bin-1).value == window.chromosome.value; ++i)
                scratch->setBin(window.orientation, --bin,
                    static_cast<std::uint16_t>(window.id.value));
            window.first_bin = BinId{bin};
            bin = window.last_bin_inclusive.value;
            for (std::uint64_t i = 0; i < config.winFlankNbins &&
                 bin+1 < config.winBinN &&
                 chromosomeForBin(bin+1).value == window.chromosome.value; ++i)
                scratch->setBin(window.orientation, ++bin,
                    static_cast<std::uint16_t>(window.id.value));
            window.last_bin_inclusive = BinId{bin};
        }

        const auto assign = [&](const SeedHitBlock& block, Orientation orientation,
                                const Seed& seed) -> bool {
            const auto bin = block.genome_start.value >> config.winBinNbits;
            const auto id = scratch->bin_to_window[static_cast<unsigned>(orientation)][bin];
            if (id == no_window) return true;
            auto& window = windows[id];
            const bool anchor = seed.multiplicity.value <= config.winAnchorMultimapNmax;
            if (!anchor && block.length < window.minimum_recordable_length) return true;
            auto& alignments = window.seeds;
            const auto start = block.read_start.value;
            const auto genome = block.genome_start.value;
            const auto junction = block.junction.value_or(JunctionId{WindowJunction::absent}).value;
            auto overlaps = [&](const WindowSeed& previous) {
                const auto previous_read = previous.read_start.value;
                return previous.mate == block.mate && previous.sjdb.encoded() == junction &&
                    genome+previous_read == previous.genome_start.value+start &&
                    ((start >= previous_read && start < previous_read+previous.length) ||
                     (start+block.length >= previous_read &&
                      start+block.length < previous_read+previous.length));
            };
            auto existing = alignments.begin();
            for (; existing != alignments.end(); ++existing)
                if (overlaps(*existing)) break;
            WindowSeed entry{block.read_start, block.genome_start, block.length,
                seed.multiplicity.value, anchor, block.mate, block.junction};
            if (existing != alignments.end()) {
                if (block.length > existing->length) {
                    // STAR assignAlignToWindow: find the insertion point while
                    // excluding the replaced entry, then shift only the interval
                    // between the two positions. erase()+insert() shifts the
                    // unaffected tail twice (even when replacing in place).
                    auto position = alignments.begin();
                    for (; position != alignments.end(); ++position)
                        if (position != existing && start < position->read_start.value)
                            break;
                    if (position > existing) --position;
                    if (position < existing)
                        std::move_backward(position, existing, existing+1);
                    else if (position > existing)
                        std::move(existing+1, position+1, existing);
                    *position = entry;
                }
                return true;
            }
            if (alignments.size() == config.seedPerWindowNmax) {
                auto minimum = read.budget.legacy_Lread+1;
                for (const auto& value : alignments)
                    if (!value.anchor) minimum = std::min(minimum, value.length);
                window.minimum_recordable_length = minimum;
                if (minimum == read.budget.legacy_Lread+1) return false;
                if (!anchor && block.length < minimum) return true;
                alignments.erase(std::remove_if(alignments.begin(), alignments.end(),
                    [minimum](const WindowSeed& value) {
                        return !value.anchor && value.length <= minimum;
                    }), alignments.end());
            }
            if (anchor || block.length > window.minimum_recordable_length) {
                auto position = alignments.begin();
                for (; position != alignments.end(); ++position)
                    if (start < position->read_start.value) break;
                alignments.insert(position, entry);
            }
            return true;
        };

        for (const auto& seed : seeds.seeds) {
            auto status = seeds.occurrences->visitBlocks(seed.occurrences,
                [&](const SeedHitBlock& block, Orientation orientation) {
                    if (!assign(block, orientation, seed)) {
                        out.markers.legacy_map_marker = 999905;
                        out.markers.normal_early_stop = true;
                        out.windows.clear();
                        return false;
                    }
                    return true;
                });
            if (std::holds_alternative<Failure>(status)) return status;
            if (out.markers.legacy_map_marker == 999905) break;
        }
        return std::monostate{};
    }
private:
    ReferenceLease reference_;
    std::string semantic_id_;
};

} // namespace

Result<std::shared_ptr<const IWindowBuilder>> makeStarWindowBuilder(
    ReferenceLease reference) {
    if (!reference) return windowError(ErrorKind::Config, "window reference is null");
    return std::shared_ptr<const IWindowBuilder>(
        std::make_shared<StarWindowBuilder>(std::move(reference)));
}

} // namespace nustar::star2711b
