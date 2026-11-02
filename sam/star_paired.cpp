#include "star_paired.hpp"

#include <algorithm>
#include <limits>

namespace nustar::star2711b {
namespace {

Failure samError(ErrorKind kind, std::string message,
                 std::optional<FragmentId> fragment = {}) {
    return Failure{kind, Stage::Project, fragment, 1, std::move(message)};
}

class SamWorkspace final : public StageWorkspace {
public:
    const WorkspaceIdentity& identity() const noexcept override { return identity_; }
    void reset() noexcept override {}
private:
    WorkspaceIdentity identity_{"star-2.7.11b-paired-sam-v1", 1, "paired"};
};

void appendTag(SamRecord& record, const char* name, char type,
               std::string value) {
    record.tags.push_back(SamAuxField{{name[0], name[1]}, type,
                                      std::move(value)});
}

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

void addCigar(std::vector<CigarOperation>& cigar, std::uint64_t length, char op) {
    if (length) cigar.push_back(CigarOperation{length, op});
}

SamRecord unmappedRecord(const PreparedReadPair& read, MateId mate,
                         const AlignmentCandidate* best, int code,
                         const ReferenceMetadata& metadata,
                         bool other_mapped, const SamConfig& config) {
    const auto index = static_cast<std::size_t>(mate);
    SamRecord record{};
    record.fragment = read.original.id;
    record.mate = mate;
    record.qname = read.output_name;
    record.flag = 0x1 | 0x4 | (index == 0 ? 0x40 : 0x80) | (read.original.legacy_read_filter == 'Y' ? 0x200 : 0);
    if (!other_mapped) record.flag |= 0x8;
    else if (best) {
        const auto strand = static_cast<unsigned>(best->orientation);
        if (strand != 1-index) record.flag |= 0x20;
    }
    record.rname = "*";
    record.rnext = "*";
    if (other_mapped && best) {
        record.rnext = metadata.chromosomes[best->chromosome.value].name;
        record.mate_position_1based = best->blocks.front().genome_start.value + 1 -
            best->chromosome_start;
    }
    record.sequence = read.original.mates[index].sequence;
    record.quality = read.quality_for_output[index];
    appendTag(record, "NH", 'i', "0");
    appendTag(record, "HI", 'i', "0");
    appendTag(record, "AS", 'i', std::to_string(best ? best->metrics.score : 0));
    appendTag(record, "nM", 'i',
              std::to_string(best ? best->metrics.mismatches : 0));
    appendTag(record, "uT", 'A', std::to_string(code));
    if(!config.read_group_by_file.empty())
        appendTag(record,"RG",'Z',config.read_group_by_file[read.original.id.file_group]);
    return record;
}

class StarPairedSamProjector final : public IPairedSamProjector {
public:
    std::string_view backendId() const noexcept override {
        return "star-2.7.11b-paired-sam-v1";
    }

    Result<std::unique_ptr<StageWorkspace>> createWorkspace(
        const PairedSamConfig& config, const CompatibilityLimits&,
        const WorkspaceBudget&) const override {
        if (config.keep_pairs || config.common.outSAMmode != "Full" ||
            std::any_of(config.common.outSAMattributes.begin(),config.common.outSAMattributes.end(),
                [](const auto& tag){return tag!="Standard" && tag!="NH" && tag!="HI" && tag!="AS" && tag!="nM" && tag!="RG" && !isSoloAttribute(tag);}))
            return samError(ErrorKind::Config,
                "paired SAM projector currently supports full text SAM with standard tags");
        return std::unique_ptr<StageWorkspace>(std::make_unique<SamWorkspace>());
    }

    Result<ProjectedOutput> projectPair(const PreparedReadPair& read,
        const CandidateSet& candidates, const SelectionResult& selection,
        const FilterDecision& decision, const ReferenceView& reference,
        const PairedSamConfig& config, const SelectionConfig& select_config,
        const PairedStitchConfig& stitch_config,
        StageWorkspace&) const override {
        ProjectedOutput out{};
        out.disposition = OutputDisposition{{false, false},
            decision.legacy_unmap_type, decision.mapping == MappingClass::Unmapped};
        if (!config.common.enabled) return out;
        const auto* best = candidates.best ?
            &candidates.candidates[candidates.best->value] : nullptr;
        if (decision.mapping == MappingClass::Unmapped) {
            if (config.common.unmapped_within)
                for (std::size_t mate = 0; mate < 2; ++mate)
                    out.records.push_back(unmappedRecord(read,
                        static_cast<MateId>(mate), best,
                        decision.legacy_unmap_type, *reference.metadata, false,config.common));
            return out;
        }
        const auto ntr = selection.ordered.size();
        const auto count = select_config.outSAMmultNmax ?
            std::min<std::size_t>(ntr, *select_config.outSAMmultNmax) : ntr;
        for (std::size_t ordinal = 0; ordinal < count; ++ordinal) {
            const auto& selected = selection.ordered[ordinal];
            const auto& c = candidates.candidates[selected.candidate.value];
            const auto gap = std::find_if(c.boundaries.begin(), c.boundaries.end(),
                [](const AlignmentBoundary& b) {
                    return b.kind == BoundaryKind::MateGap;
                });
            const bool both = gap != c.boundaries.end();
            const auto first_last = both ?
                std::size_t(gap-c.boundaries.begin()) : c.blocks.size()-1;
            const auto& chromosome = reference.metadata->chromosomes[c.chromosome.value].name;
            const auto strand = static_cast<unsigned>(c.orientation);
            const auto left_mate = strand;
            bool proper = false;
            if (both) {
                const auto& left = c.blocks.front();
                const auto& left_end = c.blocks[first_last];
                const auto& right = c.blocks[first_last+1];
                const auto& right_end = c.blocks.back();
                proper = stitch_config.protrusion_concordant_pair ||
                    (left.genome_start.value <= right.genome_start.value+
                        left.read_start.value &&
                     left_end.genome_start.value+left_end.length <=
                        right_end.genome_start.value+
                        read.mapping.view().budget.legacy_Lread-
                        right_end.read_start.value);
            }
            const auto span = c.blocks.back().genome_start.value+
                c.blocks.back().length-c.blocks.front().genome_start.value;
            const std::size_t groups = both ? 2 : 1;
            for (std::size_t group = 0; group < groups; ++group) {
                const auto begin = group == 0 ? 0 : first_last+1;
                const auto end = group == 0 ? first_last : c.blocks.size()-1;
                const auto mate = c.blocks[begin].mate;
                const auto mate_index = static_cast<std::size_t>(mate);
                SamRecord record{};
                record.fragment = read.original.id;
                record.mate = mate;
                record.alignment_ordinal = ordinal;
                record.qname = read.output_name;
                record.flag = 0x1 | (proper ? 0x2 : 0) |
                    (both ? 0 : 0x8) | (mate_index == 0 ? 0x40 : 0x80) |
                    (selected.primary ? 0 : 0x100);
                record.flag |= mate_index == 0 ? strand*0x10 : (1-strand)*0x10;
                if (both) record.flag |= mate_index == 0 ?
                    (1-strand)*0x20 : strand*0x20;
                record.rname = chromosome;
                record.position_1based = c.blocks[begin].genome_start.value+1-
                    c.chromosome_start;
                record.mapq = ntr >= 5 ? 0 : ntr >= 3 ? 1 : ntr == 2 ? 3 :
                    config.common.outSAMmapqUnique;
                const auto local_origin = c.blocks[begin].read_start.value <
                    read.mapping.view().mates[left_mate].clipped_length ? 0 :
                    read.mapping.view().mates[left_mate].clipped_length+1;
                const auto local_start = c.blocks[begin].read_start.value-local_origin;
                const auto& clipping=read.mapping.view().mates[mate_index].clipping;
                const auto trim_left=mate_index==strand ? clipping.total_5p:clipping.total_3p;
                addCigar(record.cigar, local_start+trim_left, 'S');
                for (std::size_t i = begin; i <= end; ++i) {
                    if (i > begin) {
                        const auto& before = c.blocks[i-1];
                        const auto& after = c.blocks[i];
                        const auto read_gap = after.read_start.value-
                            before.read_start.value-before.length;
                        const auto genome_gap = after.genome_start.value-
                            before.genome_start.value-before.length;
                        addCigar(record.cigar, read_gap, 'I');
                        addCigar(record.cigar, genome_gap,
                            c.boundaries[i-1].kind == BoundaryKind::Intron ? 'N' : 'D');
                    }
                    addCigar(record.cigar, c.blocks[i].length, 'M');
                }
                const auto local_end = c.blocks[end].read_start.value+
                    c.blocks[end].length-local_origin;
                addCigar(record.cigar,
                    read.original.mates[mate_index].sequence.size()-local_end-trim_left, 'S');
                if (both) {
                    record.rnext = "=";
                    record.mate_position_1based = c.blocks[group == 0 ?
                        first_last+1 : 0].genome_start.value+1-c.chromosome_start;
                    record.template_length = group == 0 ? span : -std::int64_t(span);
                } else record.rnext = "*";
                const auto& raw = read.original.mates[mate_index];
                record.sequence = mate_index == strand ? raw.sequence :
                    reverseComplement(raw.sequence);
                record.quality = mate_index == strand ?
                    read.quality_for_output[mate_index] :
                    std::string(read.quality_for_output[mate_index].rbegin(),
                                read.quality_for_output[mate_index].rend());
                record.tags.reserve(4);
                if(read.original.legacy_read_filter=='Y') record.flag|=0x200;
                appendTag(record, "NH", 'i', std::to_string(ntr));
                appendTag(record, "HI", 'i',
                    std::to_string(ordinal+config.common.outSAMattrIHstart));
                appendTag(record, "AS", 'i', std::to_string(c.metrics.score));
                appendTag(record, "nM", 'i',
                    std::to_string(c.metrics.mismatches));
                if(!config.common.read_group_by_file.empty() &&
                   std::find(config.common.outSAMattributes.begin(),config.common.outSAMattributes.end(),"RG")!=config.common.outSAMattributes.end())
                    appendTag(record,"RG",'Z',config.common.read_group_by_file[read.original.id.file_group]);
                if(!config.common.outSAMattributes.empty() && config.common.outSAMattributes!=std::vector<std::string>{"Standard"}) {
                    std::vector<SamAuxField> ordered;
                    for(const auto& tag:config.common.outSAMattributes)
                        for(auto& field:record.tags)
                            if(field.tag[0]==tag[0] && field.tag[1]==tag[1]) {ordered.push_back(std::move(field));break;}
                    record.tags=std::move(ordered);
                }
                out.records.push_back(std::move(record));
            }
        }
        if (best) {
            for (const auto& block : best->blocks)
                out.disposition.mate_mapped[static_cast<std::size_t>(block.mate)] = true;
        }
        if (!out.disposition.mate_mapped[0] ||
            !out.disposition.mate_mapped[1]) {
            out.disposition.legacy_unmap_type_after_output = 4;
            out.disposition.count_as_unmapped_all = true;
            if (config.common.unmapped_within)
                for (std::size_t mate = 0; mate < 2; ++mate)
                    if (!out.disposition.mate_mapped[mate])
                        out.records.push_back(unmappedRecord(read,
                            static_cast<MateId>(mate), best, 4,
                            *reference.metadata, true,config.common));
        }
        return out;
    }
};

} // namespace

std::shared_ptr<const IPairedSamProjector> makeStarPairedSamProjector() {
    return std::make_shared<StarPairedSamProjector>();
}

} // namespace nustar::star2711b
