#include "star_single.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <utility>

namespace nustar::star2711b {
namespace {

Failure samError(ErrorKind kind, std::string message,
                 std::optional<FragmentId> fragment = {}) {
    return Failure{kind, Stage::Project, fragment, 1, std::move(message)};
}

class SamWorkspace final : public StageWorkspace {
public:
    explicit SamWorkspace(std::span<const std::string> tags)
        : attributes(tags.begin(),tags.end()),
          emit_motif(std::find(tags.begin(),tags.end(),"jM")!=tags.end()),
          emit_intron(std::find(tags.begin(),tags.end(),"jI")!=tags.end()) {}
    const std::vector<std::string> attributes;
    const bool emit_motif, emit_intron;
    const WorkspaceIdentity& identity() const noexcept override { return identity_; }
    void reset() noexcept override {}
private:
    WorkspaceIdentity identity_{"star-2.7.11b-single-sam-v1", 1, "single"};
};

char complement(char base) {
    switch (base) {
    case 'A': return 'T'; case 'a': return 't';
    case 'C': return 'G'; case 'c': return 'g';
    case 'G': return 'C'; case 'g': return 'c';
    case 'T': return 'A'; case 't': return 'a';
    default: return base;
    }
}

std::string reverseComplement(std::string sequence) {
    std::reverse(sequence.begin(), sequence.end());
    for (auto& base : sequence) base = complement(base);
    return sequence;
}

std::span<const std::string> tagsFor(const SamConfig& config) {
    static const std::array<std::string, 4> standard{"NH", "HI", "AS", "nM"};
    if (config.outSAMattributes.empty() ||
        (config.outSAMattributes.size() == 1 && config.outSAMattributes[0] == "Standard"))
        return standard;
    return config.outSAMattributes;
}

void appendTag(SamRecord& record, const std::string& name,
               char type, std::string value) {
    record.tags.push_back(SamAuxField{{name[0], name[1]}, type,
                                      std::move(value)});
}

char baseLetter(BaseCode base) {
    static constexpr char letters[] = "ACGTN";
    const auto value = static_cast<unsigned>(base);
    return value < 5 ? letters[value] : 'N';
}

Result<std::pair<std::uint64_t, std::string>> nmAndMd(
    const PreparedSingleRead& read, const AlignmentCandidate& candidate,
    const ReferenceView& reference) {
    const auto& oriented = read.mapping.view().encoded[
        candidate.original_read_orientation == Orientation::Forward ? 0 : 2];
    std::uint64_t nm = 0, match_count = 0;
    std::string md;
    for (std::size_t i = 0; i < candidate.blocks.size(); ++i) {
        const auto& block = candidate.blocks[i];
        auto pinned = reference.sequence->pin(GenomeRange{
            block.genome_start, GenomeOffset{block.genome_start.value+block.length}});
        if (std::holds_alternative<Failure>(pinned))
            return std::get<Failure>(std::move(pinned));
        const auto& bases = std::get<SequenceWindow>(pinned).bases;
        for (std::uint64_t j = 0; j < block.length; ++j) {
            const auto read_base = oriented[block.read_start.value+j];
            const auto genome_base = bases[j];
            if (read_base != genome_base || read_base == BaseCode::N ||
                genome_base == BaseCode::N) {
                ++nm;
                md += std::to_string(match_count);
                md += baseLetter(genome_base);
                match_count = 0;
            } else ++match_count;
        }
        if (i+1 < candidate.blocks.size()) {
            const auto& next = candidate.blocks[i+1];
            const auto& boundary = candidate.boundaries[i];
            const auto read_gap = next.read_start.value-
                (block.read_start.value+block.length);
            if (boundary.kind == BoundaryKind::Insertion) nm += read_gap;
            if (boundary.kind == BoundaryKind::Deletion) {
                const auto deleted = next.genome_start.value-
                    (block.genome_start.value+block.length);
                auto deletion = reference.sequence->pin(GenomeRange{
                    GenomeOffset{block.genome_start.value+block.length},
                    GenomeOffset{next.genome_start.value}});
                if (std::holds_alternative<Failure>(deletion))
                    return std::get<Failure>(std::move(deletion));
                nm += deleted;
                md += std::to_string(match_count);
                md += '^';
                for (auto base : std::get<SequenceWindow>(deletion).bases)
                    md += baseLetter(base);
                match_count = 0;
            }
        }
    }
    md += std::to_string(match_count);
    return std::pair{nm, std::move(md)};
}

class StarSingleSamProjector final : public ISingleSamProjector {
public:
    std::string_view backendId() const noexcept override {
        return "star-2.7.11b-single-sam-v1";
    }

    Result<std::unique_ptr<StageWorkspace>> createWorkspace(
        const SamConfig& config, const CompatibilityLimits&,
        const WorkspaceBudget&) const override {
        if (config.outSAMmode != "Full" && config.outSAMmode != "NoQS" &&
            config.outSAMmode != "None")
            return samError(ErrorKind::Config, "unsupported SAM mode");
        for (const auto& tag : tagsFor(config))
            if (tag != "NH" && tag != "HI" && tag != "AS" && tag != "nM" &&
                tag != "jM" && tag != "jI" && tag != "XS" &&
                tag != "NM" && tag != "MD" && tag != "RG" && !isSoloAttribute(tag))
                return samError(ErrorKind::Config,
                                "unsupported single-end SAM attribute: " + tag);
        return std::unique_ptr<StageWorkspace>(std::make_unique<SamWorkspace>(tagsFor(config)));
    }

    Result<ProjectedOutput> projectSingle(
        const PreparedSingleRead& read, const CandidateSet& candidates,
        const SelectionResult& selected, const FilterDecision& decision,
        const ReferenceView& reference, const SamConfig& config,
        const SelectionConfig& select_config, const SingleStitchConfig&,
        StageWorkspace& workspace) const override {
        const auto& scratch=static_cast<const SamWorkspace&>(workspace);
        const auto& attributes=scratch.attributes;
        ProjectedOutput output{};
        const bool mapped = decision.mapping != MappingClass::Unmapped;
        output.disposition = OutputDisposition{{mapped, false},
            decision.legacy_unmap_type, !mapped};
        if (!config.enabled || config.outSAMmode == "None") return output;
        const auto& sequence = read.original.mate.sequence;
        const std::string_view quality = read.original.format == FastxFormat::Fastq &&
                             config.outSAMmode != "NoQS" ?
            std::string_view(read.quality_for_output) : std::string_view("*");
        if (!mapped) {
            if (!config.unmapped_within) return output;
            SamRecord record{};
            record.fragment = read.original.id;
            record.mate = MateId::First;
            record.qname = read.output_name;
            record.flag = 0x4 | (read.original.legacy_read_filter == 'Y' ? 0x200 : 0);
            record.rname = record.rnext = "*";
            record.sequence = sequence;
            record.quality = read.original.format == FastxFormat::Fastq ?
                read.quality_for_output : "*";
            const auto* best = candidates.best ?
                &candidates.candidates[candidates.best->value] : nullptr;
            appendTag(record, "NH", 'i', "0");
            appendTag(record, "HI", 'i', "0");
            appendTag(record, "AS", 'i',
                      std::to_string(best ? best->metrics.score : 0));
            appendTag(record, "nM", 'i',
                      std::to_string(best ? best->metrics.mismatches : 0));
            appendTag(record, "uT", 'A',
                      std::to_string(decision.legacy_unmap_type));
            if (!config.read_group_by_file.empty() &&
                read.original.id.file_group < config.read_group_by_file.size())
                appendTag(record, "RG", 'Z',
                          config.read_group_by_file[read.original.id.file_group]);
            output.records.push_back(std::move(record));
            return output;
        }

        const auto count = selected.ordered.size();
        const auto write_count = select_config.outSAMmultNmax ?
            std::min<std::size_t>(count, *select_config.outSAMmultNmax) : count;
        for (std::size_t i = 0; i < write_count; ++i) {
            const auto& selection = selected.ordered[i];
            const auto& candidate = candidates.candidates[selection.candidate.value];
            SamRecord record{};
            record.fragment = read.original.id;
            record.mate = MateId::First;
            record.alignment_ordinal = i;
            record.qname = read.output_name;
            record.flag = (candidate.orientation == Orientation::Reverse ? 0x10 : 0) |
                          (selection.primary ? 0 : 0x100) |
                          (read.original.legacy_read_filter == 'Y' ? 0x200 : 0);
            record.flag = (record.flag & config.outSAMflagAND) | config.outSAMflagOR;
            record.rname = reference.metadata->chromosomes[candidate.chromosome.value].name;
            record.position_1based = candidate.blocks.front().genome_start.value+1-
                candidate.chromosome_start;
            record.mapq = count >= 5 ? 0 : count >= 3 ? 1 :
                           count == 2 ? 3 : config.outSAMmapqUnique;
            record.rnext = "*";
            record.sequence = candidate.orientation == Orientation::Forward ?
                sequence : reverseComplement(sequence);
            if (candidate.orientation == Orientation::Forward || quality == "*")
                record.quality.assign(quality);
            else record.quality.assign(quality.rbegin(), quality.rend());
            const auto& mate = read.mapping.view().mates.front();
            const auto trim_left = candidate.orientation == Orientation::Forward ?
                mate.clipping.total_5p : mate.clipping.total_3p;
            const auto first = candidate.blocks.front().read_start.value;
            const auto& last = candidate.blocks.back();
            if (first+trim_left)
                record.cigar.push_back(CigarOperation{first+trim_left, 'S'});
            std::string motif = "c", intron = "i";
            bool has_junction = false;
            for (std::size_t block = 0; block < candidate.blocks.size(); ++block) {
                const auto& current = candidate.blocks[block];
                if (block) {
                    const auto& previous = candidate.blocks[block-1];
                    const auto& boundary = candidate.boundaries[block-1];
                    const auto read_gap = current.read_start.value-
                        previous.read_start.value-previous.length;
                    const auto genome_gap = current.genome_start.value-
                        previous.genome_start.value-previous.length;
                    if (read_gap) record.cigar.push_back(CigarOperation{read_gap, 'I'});
                    if (boundary.kind == BoundaryKind::Intron || boundary.annotation) {
                        record.cigar.push_back(CigarOperation{genome_gap, 'N'});
                        if (scratch.emit_motif) {
                            motif += ',';
                            motif += std::to_string(static_cast<int>(boundary.motif)+
                                (boundary.annotation ? 20 : 0));
                        }
                        if (scratch.emit_intron) {
                            intron += ',';
                            intron += std::to_string(previous.genome_start.value+
                                previous.length+1-candidate.chromosome_start);
                            intron += ',';
                            intron += std::to_string(current.genome_start.value-
                                candidate.chromosome_start);
                        }
                        has_junction = true;
                    } else if (genome_gap)
                        record.cigar.push_back(CigarOperation{genome_gap, 'D'});
                }
                record.cigar.push_back(CigarOperation{current.length, 'M'});
            }
            const auto suffix = mate.original_length-last.read_start.value-
                last.length-trim_left;
            if (suffix) record.cigar.push_back(CigarOperation{suffix, 'S'});
            if (!has_junction) { motif = "c,-1"; intron = "i,-1"; }
            std::optional<std::pair<std::uint64_t, std::string>> mismatch_tags;
            record.tags.reserve(attributes.size());
            for (const auto& tag : attributes) {
                if (tag == "NH") appendTag(record, tag, 'i', std::to_string(count));
                else if (tag == "HI") appendTag(record, tag, 'i',
                    std::to_string(i+config.outSAMattrIHstart));
                else if (tag == "AS") appendTag(record, tag, 'i',
                    std::to_string(candidate.metrics.score));
                else if (tag == "nM") appendTag(record, tag, 'i',
                    std::to_string(candidate.metrics.mismatches));
                else if (tag == "jM") appendTag(record, tag, 'B', motif);
                else if (tag == "jI") appendTag(record, tag, 'B', intron);
                else if (tag == "XS") {
                    if (candidate.metrics.motif_strand == Strand::Plus)
                        appendTag(record, tag, 'A', "+");
                    else if (candidate.metrics.motif_strand == Strand::Minus)
                        appendTag(record, tag, 'A', "-");
                } else if (tag == "NM" || tag == "MD") {
                    if (!mismatch_tags) {
                        auto calculated = nmAndMd(read, candidate, reference);
                        if (std::holds_alternative<Failure>(calculated))
                            return std::get<Failure>(std::move(calculated));
                        mismatch_tags = std::get<std::pair<std::uint64_t, std::string>>(
                            std::move(calculated));
                    }
                    appendTag(record, tag, tag == "NM" ? 'i' : 'Z',
                        tag == "NM" ? std::to_string(mismatch_tags->first) :
                                      mismatch_tags->second);
                } else if (tag == "RG" &&
                           read.original.id.file_group < config.read_group_by_file.size())
                    appendTag(record, tag, 'Z',
                        config.read_group_by_file[read.original.id.file_group]);
            }
            output.records.push_back(std::move(record));
        }
        return output;
    }
};

} // namespace

std::shared_ptr<const ISingleSamProjector> makeStarSingleSamProjector() {
    return std::make_shared<StarSingleSamProjector>();
}

} // namespace nustar::star2711b
