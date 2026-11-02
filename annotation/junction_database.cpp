#include "junction_database.hpp"
#include <algorithm>
#include <tuple>

namespace nustar::star2711b {

namespace {
bool validBase(BaseCode base) { return base <= BaseCode::T; }

Motif motifOf(BaseCode a, BaseCode b, BaseCode c, BaseCode d) {
    const auto A=BaseCode::A, C=BaseCode::C, G=BaseCode::G, T=BaseCode::T;
    if (a==G && b==T && c==A && d==G) return Motif::GT_AG;
    if (a==C && b==T && c==A && d==C) return Motif::CT_AC;
    if (a==G && b==C && c==A && d==G) return Motif::GC_AG;
    if (a==C && b==T && c==G && d==C) return Motif::CT_GC;
    if (a==A && b==T && c==A && d==C) return Motif::AT_AC;
    if (a==G && b==T && c==A && d==T) return Motif::GT_AT;
    return Motif::Noncanonical;
}

struct PreparedProposal {
    std::size_t source;
    GenomeRange normalized;
    GenomeRange second_pass;
    GenomeRange reported;
    Motif motif;
    Strand strand;
    std::uint8_t left, right, priority;
};

bool preferred(const PreparedProposal& candidate, const PreparedProposal& current) {
    if (candidate.priority != current.priority) return candidate.priority > current.priority;
    if ((candidate.strand != Strand::Unknown) != (current.strand != Strand::Unknown))
        return candidate.strand != Strand::Unknown;
    if ((candidate.motif != Motif::Noncanonical) != (current.motif != Motif::Noncanonical))
        return candidate.motif != Motif::Noncanonical;
    if (candidate.left != current.left) return candidate.left < current.left;
    return candidate.source < current.source;
}

} // namespace

Result<PreparedJunctionDatabase> GenomeJunctionPreparer::prepare(
    std::span<const JunctionProposal> proposals, const ReferenceView& reference,
    const AnnotationConfig& config) {
    PreparedJunctionDatabase result{};
    result.overhang = config.sjdbOverhang;
    if (proposals.size() > config.limitSjdbInsertNsj)
        return Failure{ErrorKind::Capacity, Stage::Sjdb, std::nullopt, 1,
                       "junction count exceeds limitSjdbInsertNsj"};
    std::vector<PreparedProposal> work;
    work.reserve(proposals.size());
    for (std::size_t i = 0; i < proposals.size(); ++i) {
        const auto& proposal = proposals[i];
        if (proposal.chromosome.value >= reference.metadata->chromosomes.size())
            return Failure{ErrorKind::Input, Stage::Sjdb, std::nullopt, 1, "invalid junction chromosome"};
        const auto& chr = reference.metadata->chromosomes[proposal.chromosome.value];
        const auto begin = proposal.intron.begin.value, end = proposal.intron.end.value;
        if (begin < chr.start.value || end > chr.start.value + chr.length || end <= begin)
            return Failure{ErrorKind::Input, Stage::Sjdb, std::nullopt, 1,
                "junction is outside chromosome: " + chr.name + " " +
                std::to_string(begin - chr.start.value + 1) + "-" +
                std::to_string(end - chr.start.value)};
        auto pinned = reference.sequence->pin({chr.start, GenomeOffset{chr.start.value + chr.length}});
        if (std::holds_alternative<Failure>(pinned)) return std::get<Failure>(std::move(pinned));
        const auto window = std::get<SequenceWindow>(std::move(pinned));
        const auto bases = window.bases;
        const auto local_begin = begin - chr.start.value, local_end = end - chr.start.value;
        const auto motif = motifOf(bases[local_begin], bases[local_begin+1],
                                   bases[local_end-2], bases[local_end-1]);
        std::uint8_t left = 0, right = 0;
        while (left < 255 && local_begin > left && local_end > left &&
               validBase(bases[local_begin-1-left]) &&
               bases[local_begin-1-left] == bases[local_end-1-left]) ++left;
        while (right < 255 && local_begin+right < chr.length &&
               local_end+right < chr.length && validBase(bases[local_begin+right]) &&
               bases[local_begin+right] == bases[local_end+right]) ++right;
        const GenomeRange normalized{GenomeOffset{begin-left}, GenomeOffset{end-left}};
        const GenomeRange reported{GenomeOffset{begin}, GenomeOffset{end}};
        work.push_back(PreparedProposal{i, normalized,
            motif == Motif::Noncanonical ? normalized : reported,
            reported, motif, proposal.strand, left, right, proposal.legacy_priority});
    }
    std::stable_sort(work.begin(), work.end(), [](const auto& a, const auto& b) {
        return std::tie(a.normalized.begin.value, a.normalized.end.value, a.strand) <
               std::tie(b.normalized.begin.value, b.normalized.end.value, b.strand);
    });
    std::vector<PreparedProposal> first_pass;
    for (const auto& item : work) {
        if (!first_pass.empty() &&
            first_pass.back().normalized.begin.value == item.normalized.begin.value &&
            first_pass.back().normalized.end.value == item.normalized.end.value &&
            first_pass.back().strand == item.strand) {
            if (preferred(item, first_pass.back())) first_pass.back() = item;
        } else first_pass.push_back(item);
    }
    std::stable_sort(first_pass.begin(), first_pass.end(), [](const auto& a, const auto& b) {
        return std::tie(a.second_pass.begin.value, a.second_pass.end.value) <
               std::tie(b.second_pass.begin.value, b.second_pass.end.value);
    });
    std::vector<PreparedProposal> final;
    for (const auto& item : first_pass) {
        if (!final.empty() &&
            final.back().second_pass.begin.value == item.second_pass.begin.value &&
            final.back().second_pass.end.value == item.second_pass.end.value) {
            if (preferred(item, final.back())) final.back() = item;
        } else final.push_back(item);
    }
    for (const auto& item : final) {
        const auto donor = item.reported.begin.value;
        const auto acceptor = item.reported.end.value;
        // STAR copies flanks from the padded genome, including chromosome padding.
        Strand strand = item.strand;
        if (strand == Strand::Unknown && item.motif != Motif::Noncanonical)
            strand = static_cast<int>(item.motif) % 2 ? Strand::Plus : Strand::Minus;
        const auto id = result.entries.size();
        result.entries.push_back(SjdbEntry{JunctionId{id}, item.second_pass,
            item.reported, strand, item.motif, item.left, item.right,
            GenomeOffset{donor - config.sjdbOverhang}, GenomeOffset{acceptor}});
        result.winning_proposal_indices.push_back(item.source);
    }
    return result;
}

} // namespace nustar::star2711b
