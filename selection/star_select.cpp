#include "star_select.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace nustar::star2711b {
namespace {

Failure selectionError(ErrorKind kind, std::string message,
                       std::optional<FragmentId> fragment = {}) {
    return Failure{kind, Stage::Select, fragment, 1, std::move(message)};
}

class StarAlignmentSelector final : public IAlignmentSelector {
public:
    Status validateConfig(const SelectionConfig& config) const override {
        if (config.outFilterMultimapNmax == 0 ||
            (config.outMultimapperOrder != "Old_2.4" &&
             config.outMultimapperOrder != "Random") ||
            (config.outSAMprimaryFlag != "OneBestScore" &&
             config.outSAMprimaryFlag != "AllBestScore"))
            return selectionError(ErrorKind::Config, "unsupported STAR selection state");
        return std::monostate{};
    }

    Status select(const MappingRead& read, const CandidateSet& candidates,
                  const SelectionConfig& config, LegacyRandomState& random,
                  SelectionResult& out) const override {
        out = SelectionResult{};
        if (candidates.windows.empty() || !candidates.best) return std::monostate{};
        const auto best_score = candidates.candidates[candidates.best->value].metrics.score;
        out.best = candidates.best;
        out.best_score = best_score;
        std::vector<CandidateId> ordered;
        for (const auto& window : candidates.windows)
            for (const auto id : window.ordered) {
                const auto& candidate = candidates.candidates[id.value];
                if (std::int64_t(candidate.metrics.score)+
                    config.outFilterMultimapScoreRange >= best_score)
                    ordered.push_back(id);
            }
        out.candidate_count = ordered.size();
        if (ordered.empty() || ordered.size() > config.outFilterMultimapNmax)
            return std::monostate{};

        std::size_t best_count = 0;
        const bool random_order = config.outMultimapperOrder == "Random";
        if (ordered.size() > 1 && (random_order || config.outSAMmultNmax)) {
            for (std::size_t i = 0; i < ordered.size(); ++i)
                if (candidates.candidates[ordered[i].value].metrics.score == best_score)
                    std::swap(ordered[i], ordered[best_count++]);
        }
        if (random_order) {
            const auto shuffle_group = [&](std::size_t start, std::size_t count) {
                for (std::size_t i = count; i-- > 1;) {
                    const auto picked = static_cast<std::size_t>(
                        random.uniform_0_1(random.engine)*i+0.5);
                    std::swap(ordered[start+i], ordered[start+picked]);
                }
            };
            shuffle_group(0, best_count);
            shuffle_group(best_count, ordered.size()-best_count);
        }
        out.normalized_coordinates_valid = true;
        out.ordered.reserve(ordered.size());
        for (std::size_t i = 0; i < ordered.size(); ++i) {
            const auto id = ordered[i];
            const auto& candidate = candidates.candidates[id.value];
            const bool primary = ordered.size() == 1 ||
                (config.outSAMprimaryFlag == "AllBestScore" ?
                    candidate.metrics.score == best_score :
                    (random_order || config.outSAMmultNmax ? i == 0 :
                        id.value == candidates.best->value));
            const auto original_start = candidate.original_read_orientation ==
                Orientation::Forward ? candidate.read_start.value :
                read.budget.legacy_Lread-candidate.read_start.value-
                    candidate.read_length;
            out.ordered.push_back(SelectedAlignment{id, primary,
                ReadOffset{original_start}, candidate.genome_start.value-
                    candidate.chromosome_start});
        }
        return std::monostate{};
    }
};

} // namespace

std::shared_ptr<const IAlignmentSelector> makeStarAlignmentSelector() {
    return std::make_shared<StarAlignmentSelector>();
}

Result<FilterDecision> filterRead(const MappingRead& read,
                                  const CandidateSet& candidates,
                                  const SelectionResult& selection,
                                  const ReadFilterConfig& config) {
    const auto unmapped = [](UnmappedReason reason, int code) {
        return FilterDecision{MappingClass::Unmapped, reason, code};
    };
    if (candidates.windows.empty() || !candidates.best)
        return unmapped(UnmappedReason::NoGoodWindow, 0);
    const auto& best = candidates.candidates[candidates.best->value];
    const auto score_ratio_threshold = static_cast<Score>(
        config.outFilterScoreMinOverLread*(read.budget.legacy_Lread-1));
    const auto match_ratio_threshold = static_cast<std::uint64_t>(
        config.outFilterMatchNminOverLread*(read.budget.legacy_Lread-1));
    if (best.metrics.score < config.outFilterScoreMin ||
        best.metrics.score < score_ratio_threshold ||
        best.metrics.matches < config.outFilterMatchNmin ||
        best.metrics.matches < match_ratio_threshold)
        return unmapped(UnmappedReason::TooShort, 1);
    if (best.metrics.mismatches > read.budget.mismatch_limit ||
        (best.read_length > 0 &&
         double(best.metrics.mismatches)/double(best.read_length) >
             config.outFilterMismatchNoverLmax))
        return unmapped(UnmappedReason::TooManyMismatches, 2);
    if (!selection.normalized_coordinates_valid)
        return unmapped(UnmappedReason::TooManyLoci, 3);
    return FilterDecision{selection.candidate_count > 1 ? MappingClass::Multiple :
                          MappingClass::Unique, UnmappedReason::None,
                          selection.candidate_count > 1 ? -2 : -1};
}

} // namespace nustar::star2711b
