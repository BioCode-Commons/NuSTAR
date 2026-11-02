#include "star_junction.hpp"

#include <algorithm>
#include <limits>
#include <map>

namespace nustar::star2711b {
namespace {

using Key = std::pair<std::uint64_t, std::uint32_t>;

Key keyOf(const JunctionKey& key) {
    return {key.start.value, key.gap};
}

Failure sjError(ErrorKind kind, std::string message,
                std::optional<FragmentId> fragment = {}) {
    return Failure{kind, Stage::SjReduce, fragment, 1, std::move(message)};
}

std::size_t motifGroup(Motif motif) {
    return (static_cast<unsigned>(motif)+1)/2;
}

std::uint64_t distance(std::uint64_t a, std::uint64_t b) {
    return a >= b ? a-b : b-a;
}

} // namespace

SjOutputConfig defaultStarSjOutputConfig() {
    SjOutputConfig config{};
    config.enabled = true;
    config.unique_reads_only = false;
    config.outSJfilterOverhangMin = {30, 12, 12, 12};
    config.outSJfilterCountUniqueMin = {3, 1, 1, 1};
    config.outSJfilterCountTotalMin = {3, 1, 1, 1};
    config.outSJfilterDistToOtherSJmin = {10, 0, 5, 10};
    config.outSJfilterIntronMaxVsReadN = {50000, 100000, 200000};
    config.limitOutSJoneRead = 1000;
    config.limitOutSJcollapsed = 1000000;
    return config;
}

Result<FragmentJunctions> collectJunctions(const MappingRead& read,
    const CandidateSet& candidates, const SelectionResult& selected,
    const FilterDecision& decision, const SjOutputConfig& config) {
    FragmentJunctions out{read.fragment, {}};
    if (!config.enabled || decision.mapping == MappingClass::Unmapped ||
        (config.unique_reads_only && selected.ordered.size() != 1))
        return out;
    std::map<Key, std::size_t> seen;
    for (const auto& chosen : selected.ordered) {
        if (chosen.candidate.value >= candidates.candidates.size())
            return sjError(ErrorKind::Input, "selected candidate missing for SJ",
                           read.fragment);
        const auto& candidate = candidates.candidates[chosen.candidate.value];
        for (std::size_t i = 0; i < candidate.boundaries.size(); ++i) {
            const auto& boundary = candidate.boundaries[i];
            if (boundary.kind != BoundaryKind::Intron) continue;
            if (i+1 >= candidate.blocks.size())
                return sjError(ErrorKind::Input, "SJ boundary lacks adjacent blocks",
                               read.fragment);
            const auto& left = candidate.blocks[i];
            const auto& right = candidate.blocks[i+1];
            const auto start = left.genome_start.value+left.length;
            if (right.genome_start.value < start ||
                right.genome_start.value-start >
                    std::numeric_limits<std::uint32_t>::max())
                return sjError(ErrorKind::Input, "invalid SJ gap", read.fragment);
            const auto gap = static_cast<std::uint32_t>(
                right.genome_start.value-start);
            const auto overhang = static_cast<std::uint16_t>(
                std::min(left.length, right.length));
            const auto key = Key{start, gap};
            const auto found = seen.find(key);
            if (found != seen.end()) {
                out.ordered[found->second].overhang =
                    std::max(out.ordered[found->second].overhang, overhang);
                continue;
            }
            seen.emplace(key, out.ordered.size());
            const auto motif_value = static_cast<unsigned>(boundary.motif);
            const auto strand = motif_value == 0 ? Strand::Unknown :
                static_cast<Strand>((motif_value+1)%2+1);
            out.ordered.push_back(JunctionObservation{
                JunctionKey{GenomeOffset{start}, gap}, strand,
                boundary.motif, boundary.annotation, overhang,
                selected.ordered.size() == 1});
            if (out.ordered.size() > config.limitOutSJoneRead)
                return sjError(ErrorKind::Capacity,
                    "junctions exceed limitOutSJoneRead", read.fragment);
        }
    }
    return out;
}

struct JunctionAccumulator::Impl {
    std::uint32_t worker;
    SjOutputConfig config;
    std::vector<JunctionObservation> observations;
    bool finished = false;
};

JunctionAccumulator::JunctionAccumulator(std::uint32_t worker,
    const SjOutputConfig& config)
    : impl_(std::make_unique<Impl>(Impl{worker, config, {}, false})) {}

JunctionAccumulator::~JunctionAccumulator() = default;

Status JunctionAccumulator::add(const FragmentJunctions& fragment) {
    if (impl_->finished)
        return sjError(ErrorKind::Input, "SJ accumulator already finished",
                       fragment.fragment);
    impl_->observations.insert(impl_->observations.end(),
        fragment.ordered.begin(), fragment.ordered.end());
    return std::monostate{};
}

Result<WorkerJunctions> JunctionAccumulator::finish() {
    if (impl_->finished)
        return sjError(ErrorKind::Input, "SJ accumulator already finished");
    impl_->finished = true;
    std::map<Key, JunctionAggregate> collapsed;
    for (const auto& observation : impl_->observations) {
        const auto key = keyOf(observation.key);
        auto [it, inserted] = collapsed.try_emplace(key,
            JunctionAggregate{observation.key, observation.strand,
                observation.motif, observation.annotation, 0, 0,
                observation.overhang, observation.overhang});
        auto& aggregate = it->second;
        if (!inserted && aggregate.motif != observation.motif)
            return sjError(ErrorKind::Internal,
                "different motifs for the same SJ key");
        if (!inserted && aggregate.annotation < observation.annotation)
            return sjError(ErrorKind::Internal,
                "SJ annotation status increased during collapse");
        if (observation.unique_fragment) ++aggregate.unique_count;
        else ++aggregate.multiple_count;
        aggregate.max_left_overhang = std::max(aggregate.max_left_overhang,
                                               observation.overhang);
        aggregate.max_right_overhang = std::max(aggregate.max_right_overhang,
                                                observation.overhang);
    }
    WorkerJunctions out{impl_->worker, {}};
    out.sorted.reserve(collapsed.size());
    for (auto& [key, value] : collapsed) out.sorted.push_back(value);
    return out;
}

Result<FinalJunctions> finalizeJunctions(std::span<const WorkerJunctions> workers,
    const SjOutputConfig& config) {
    FinalJunctions out{};
    if (!config.enabled) return out;
    std::map<Key, JunctionAggregate> merged;
    for (const auto& worker : workers)
        for (const auto& observation : worker.sorted) {
            const auto key = keyOf(observation.key);
            auto [it, inserted] = merged.try_emplace(key, observation);
            if (inserted) continue;
            auto& aggregate = it->second;
            if (aggregate.motif != observation.motif ||
                aggregate.annotation < observation.annotation)
                return sjError(ErrorKind::Internal,
                    "incompatible SJ annotations across workers");
            aggregate.unique_count += observation.unique_count;
            aggregate.multiple_count += observation.multiple_count;
            aggregate.max_left_overhang = std::max(aggregate.max_left_overhang,
                                                   observation.max_left_overhang);
            aggregate.max_right_overhang = std::max(aggregate.max_right_overhang,
                                                    observation.max_right_overhang);
        }
    std::vector<JunctionAggregate> supported;
    supported.reserve(merged.size());
    for (auto& [key, sj] : merged) {
        const auto motif = motifGroup(sj.motif);
        const auto total = std::uint64_t(sj.unique_count)+sj.multiple_count;
        const bool enough = sj.annotation > 0 ||
            ((sj.unique_count >= std::uint32_t(config.outSJfilterCountUniqueMin[motif]) ||
              total >= std::uint32_t(config.outSJfilterCountTotalMin[motif])) &&
             sj.max_left_overhang >= config.outSJfilterOverhangMin[motif] &&
             sj.max_right_overhang >= config.outSJfilterOverhangMin[motif] &&
             (total > config.outSJfilterIntronMaxVsReadN.size() ||
              sj.key.gap <= std::uint32_t(
                  config.outSJfilterIntronMaxVsReadN[total-1])));
        if (enough) supported.push_back(sj);
    }
    std::vector<bool> passed(supported.size(), false);
    for (std::size_t i = 0; i < supported.size(); ++i) {
        const auto& sj = supported[i];
        const auto previous = i == 0 ? 0 : supported[i-1].key.start.value;
        const auto next = i+1 == supported.size() ?
            std::numeric_limits<std::uint64_t>::max() :
            supported[i+1].key.start.value;
        const auto min_distance = std::min(distance(sj.key.start.value, previous),
                                           distance(next, sj.key.start.value));
        passed[i] = min_distance >= std::uint32_t(
            config.outSJfilterDistToOtherSJmin[motifGroup(sj.motif)]);
    }
    std::vector<std::pair<std::uint64_t, std::size_t>> acceptors;
    acceptors.reserve(supported.size());
    for (std::size_t i = 0; i < supported.size(); ++i)
        acceptors.emplace_back(supported[i].key.start.value+
                               supported[i].key.gap, i);
    std::sort(acceptors.begin(), acceptors.end(),
        [](const auto& a, const auto& b) { return a.first < b.first; });
    for (std::size_t i = 0; i < acceptors.size(); ++i) {
        const auto index = acceptors[i].second;
        const auto& sj = supported[index];
        if (sj.annotation) {
            passed[index] = true;
            continue;
        }
        const auto previous = i == 0 ? 0 : acceptors[i-1].first;
        const auto next = i+1 == acceptors.size() ?
            std::numeric_limits<std::uint64_t>::max() : acceptors[i+1].first;
        const auto min_distance = std::min(distance(acceptors[i].first, previous),
                                           distance(next, acceptors[i].first));
        passed[index] = passed[index] && min_distance >= std::uint32_t(
            config.outSJfilterDistToOtherSJmin[motifGroup(sj.motif)]);
    }
    for (std::size_t i = 0; i < supported.size(); ++i) {
        out.trace.push_back(JunctionFilterTrace{
            supported[i].key, true, passed[i], passed[i], passed[i]});
        if (passed[i]) out.accepted.push_back(supported[i]);
    }
    return out;
}

} // namespace nustar::star2711b
