#include "../../pipeline/single_end/pipeline.hpp"
#include "../../pipeline/paired_end/pipeline.hpp"
#include "../../reference/genome.hpp"
#include "../../reads/fastq_prepare.hpp"
#include "../../runtime/defaults.hpp"
#include "../../sam/star_single.hpp"
#include "../../sam/star_paired.hpp"
#include "../../selection/star_select.hpp"
#include "../../splice_junction/star_junction.hpp"
#include "../../statistics/interface.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace nustar::star2711b;

namespace {
void expect(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
template<class T> T take(Result<T> value) {
    if (const auto* error = std::get_if<Failure>(&value))
        throw std::runtime_error(error->message);
    return std::get<T>(std::move(value));
}
// No mapping is needed to reject an invalid configuration, even for empty input.
class UnusedMapper final : public ISingleCandidateMapper, public IPairedCandidateMapper {
public:
    std::string_view backendId() const noexcept override { return "unused"; }
    Result<std::unique_ptr<StageWorkspace>> createWorkspace(const WorkspaceBudget&) const override {
        throw std::runtime_error("mapping workspace created during configuration validation");
    }
    Status mapSingle(const SingleMappingRead&, StageWorkspace&, CandidateSet&) const override {
        throw std::runtime_error("single mapping called during configuration validation");
    }
    Status mapPair(const PairedMappingRead&, StageWorkspace&, CandidateSet&) const override {
        throw std::runtime_error("paired mapping called during configuration validation");
    }
};
// Existing alternative selectors need not implement the new optional hook.
class AlternativeSelector final : public IAlignmentSelector {
public:
    Status select(const MappingRead&, const CandidateSet&, const SelectionConfig&,
                  LegacyRandomState&, SelectionResult& out) const override {
        out = {};
        return std::monostate{};
    }
};
}

int main() {
    try {
        auto bases = std::make_shared<const std::vector<BaseCode>>(256, BaseCode::A);
        auto storage = std::make_shared<GenomeStorage>();
        storage->sequence_storage = bases;
        storage->forward = *bases;
        auto reference = std::make_shared<GenomeReferenceData>(storage,
            std::vector<SjdbEntry>{}, std::vector<BaseCode>{}, 0);
        auto mapper = std::make_shared<UnusedMapper>();
        auto selector = makeStarAlignmentSelector();
        const SingleKernelBundle single{mapper, selector, makeStarSingleSamProjector()};
        const PairedKernelBundle paired{mapper, {}, {}, selector, makeStarPairedSamProjector()};
        auto validatePipelines = [&](const SelectionConfig& selection, bool accepted) {
            SinglePipelineConfig se{};
            PairedPipelineConfig pe{};
            se.select = pe.select = selection;
            const auto a = composeSinglePipeline(reference, single, se);
            const auto b = composePairedPipeline(reference, paired, pe);
            expect(std::holds_alternative<Failure>(a) != accepted,
                   "SE configuration acceptance differs");
            expect(std::holds_alternative<Failure>(b) != accepted,
                   "PE configuration acceptance differs");
            if (!accepted) {
                expect(std::get<Failure>(a).kind == ErrorKind::Config &&
                       std::get<Failure>(b).kind == ErrorKind::Config,
                       "configuration failure kind differs");
            }
        };
        for (const auto* order : {"Old_2.4", "Random"})
            for (const auto* primary : {"OneBestScore", "AllBestScore"}) {
                auto config = defaultStarSelectionConfig();
                config.outMultimapperOrder = order;
                config.outSAMprimaryFlag = primary;
                validatePipelines(config, true);
            }
        auto config = defaultStarSelectionConfig();
        config.outFilterMultimapNmax = 0;
        validatePipelines(config, false);
        config = defaultStarSelectionConfig();
        config.outMultimapperOrder = "unsupported";
        validatePipelines(config, false);
        config = defaultStarSelectionConfig();
        config.outSAMprimaryFlag = "unsupported";
        validatePipelines(config, false);
        expect(!std::holds_alternative<Failure>(AlternativeSelector{}.validateConfig(config)),
               "default validation hook rejected alternative selector configuration");

        config = defaultStarSelectionConfig();
        take(selector->validateConfig(config));
        auto prepared = prepareSingleFastq("@splice", "AAAAAA", "IIIIII", 0);
        const auto& read = prepared.mapping.view();
        CandidateSet candidates{};
        AlignmentCandidate candidate{};
        candidate.id = CandidateId{0};
        candidate.window = WindowId{0};
        candidate.read_length = 6;
        candidate.genome_start = GenomeOffset{10};
        candidate.blocks = {
            {ReadOffset{0}, GenomeOffset{10}, 3, MateId::First, {}},
            {ReadOffset{3}, GenomeOffset{33}, 3, MateId::First, {}}};
        candidate.boundaries = {{BoundaryKind::Intron, Motif::GT_AG, 1, Strand::Plus, 0, 0}};
        candidate.metrics.score = 6;
        candidate.metrics.matches = 6;
        candidates.candidates.push_back(candidate);
        candidates.windows.push_back({WindowId{0}, {CandidateId{0}}});
        candidates.best = CandidateId{0};
        SelectionResult selected{};
        LegacyRandomState random{};
        take(selector->select(read, candidates, config, random, selected));
        const auto decision = take(filterRead(read, candidates, selected, defaultStarReadFilterConfig()));
        expect(decision.mapping == MappingClass::Unique, "valid splice was not uniquely mapped");
        const auto sj_config = defaultStarSjOutputConfig();
        const auto junctions = take(collectJunctions(read, candidates, selected, decision, sj_config));
        expect(junctions.ordered.size() == 1 && junctions.ordered[0].key.start.value == 13 &&
               junctions.ordered[0].key.gap == 20, "adjacent splice blocks produced wrong junction");
        const auto statistics = take(deriveStatistics(read, candidates, selected, decision, {}));
        expect(statistics.mapped_unique == 1 && statistics.mapped_bases == 6 &&
               statistics.mapped_portion == 1.0 && statistics.annotated_splices == 1,
               "unique read statistics changed");
        JunctionAccumulator accumulator(0, sj_config);
        take(accumulator.add(junctions));
        auto next_fragment = junctions;
        next_fragment.fragment.ordinal = 1;
        take(accumulator.add(next_fragment));
        const auto collapsed = take(accumulator.finish());
        expect(collapsed.sorted.size() == 1 && collapsed.sorted[0].unique_count == 2,
               "accumulator lost committed junction observations");

        for (auto empty : {prepareSingleFastq("@empty", "", "", 1),
                           prepareSingleFastq("@clipped", "AAAA", "IIII", 2, {4, 0})}) {
            CandidateSet none{};
            take(selector->select(empty.mapping.view(), none, config, random, selected));
            const auto unmapped = take(filterRead(empty.mapping.view(), none, selected,
                                                  defaultStarReadFilterConfig()));
            expect(unmapped.mapping == MappingClass::Unmapped, "empty read mapped");
            expect(take(collectJunctions(empty.mapping.view(), none, selected, unmapped,
                                        sj_config)).ordered.empty(), "empty read emitted SJ");
            const auto stats = take(deriveStatistics(empty.mapping.view(), none, selected, unmapped, {}));
            expect(stats.mapped_unique == 0 && std::isfinite(stats.mapped_portion),
                   "empty read entered unique statistics division");
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
