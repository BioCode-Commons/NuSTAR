#include "../../sam/star_single.hpp"
#include "../../reads/fastq_prepare.hpp"
#include <array>
#include <iostream>
#include <stdexcept>
#include <utility>
using namespace nustar::star2711b;
namespace {
template<class T> T take(Result<T> result) {
    if (auto* failure=std::get_if<Failure>(&result)) throw std::runtime_error(failure->message);
    return std::get<T>(std::move(result));
}
void require(bool condition,const char* message) { if(!condition) throw std::runtime_error(message); }
class Sequence final : public IReferenceSequence {
    std::array<BaseCode,100> bases{};
public:
    Result<SequenceWindow> pin(GenomeRange range) const override {
        return SequenceWindow{range,std::span(bases).subspan(range.begin.value,
            range.end.value-range.begin.value),{}};
    }
};
std::string tag(const SamRecord& record,std::string_view name) {
    for(const auto& field:record.tags)
        if(std::string_view(field.tag.data(),2)==name) return field.encoded_value;
    throw std::runtime_error("missing tag: "+std::string(name));
}
}
int main() {
    try {
        auto read=prepareSingleFastq("@probe", "AAAAAA", "ABCDEF",0);
        require(read.original.mate.sequence=="AAAAAA" && *read.original.mate.quality=="ABCDEF"
            && read.quality_for_output=="ABCDEF","prepared read lost owned strings");
        auto pair=preparePairedFastq("@pair/1",std::string(150,'A'),std::string(150,'I'),
            "@pair/2",std::string(150,'C'),std::string(150,'J'),1);
        require(pair.original.mates[1].sequence==std::string(150,'C') &&
            *pair.original.mates[0].quality==std::string(150,'I') &&
            pair.quality_for_output[1]==std::string(150,'J'),"paired preparation lost strings");
        Sequence sequence;
        ReferenceMetadata metadata{{Chromosome{ChromosomeId{0},"chr1",GenomeOffset{0},100}},GenomeOffset{100},"test"};
        ReferenceView reference{ReferenceGeneration{1},"test",&metadata,&sequence,{}};
        CandidateSet candidates{};
        AlignmentCandidate candidate{};
        candidate.chromosome=ChromosomeId{0};
        candidate.orientation=candidate.original_read_orientation=Orientation::Forward;
        candidate.blocks={{{0},{0},3,MateId::First,{}},{{3},{7},3,MateId::First,{}}};
        candidate.boundaries={{BoundaryKind::Intron,Motif::GT_AG,1,Strand::Plus,0,0}};
        candidate.metrics.score=6;
        candidate.metrics.motif_strand=Strand::Plus;
        candidates.candidates.push_back(candidate);candidates.best=CandidateId{0};
        SelectionResult selected{};
        selected.ordered.push_back({CandidateId{0},true,ReadOffset{0},0});
        auto projector=makeStarSingleSamProjector();
        SamConfig config{};config.enabled=true;config.outSAMmode="Full";
        config.outSAMattrIHstart=1;config.outSAMmapqUnique=255;config.unmapped_within=true;
        auto project=[&](FilterDecision decision) {
            auto workspace=take(projector->createWorkspace(config,{},{}));
            return take(projector->projectSingle(read,candidates,selected,decision,
                reference,config,{}, {},*workspace));
        };
        const FilterDecision mapped{MappingClass::Unique,UnmappedReason::None,-1};
        for(const auto& attributes:std::vector<std::vector<std::string>>{
                {"Standard"},{"jM"},{"jI"},{"MD","jI","jM","NM","XS","NH"}}) {
            config.outSAMattributes=attributes;
            auto output=project(mapped);const auto& record=output.records.at(0);
            require(record.sequence=="AAAAAA" && record.quality=="ABCDEF","mapped strings differ");
            require(record.cigar.size()==3 && record.cigar[1].operation=='N' &&
                record.cigar[1].length==4,"spliced CIGAR differs");
            if(attributes[0]=="Standard") require(record.tags.size()==4,"standard tag list changed");
            else for(std::size_t i=0;i<attributes.size();++i) {
                require(std::string(record.tags[i].tag.data(),2)==attributes[i],"tag order changed");
                if(attributes[i]=="jM") require(tag(record,"jM")=="c,21","annotated motif differs");
                if(attributes[i]=="jI") require(tag(record,"jI")=="i,4,7","junction coordinates differ");
                if(attributes[i]=="NM") require(tag(record,"NM")=="0","NM differs");
                if(attributes[i]=="MD") require(tag(record,"MD")=="6","MD differs");
            }
        }
        config.outSAMattributes={"jM","jI"};
        candidates.candidates[0].blocks.resize(1);candidates.candidates[0].blocks[0].length=6;
        candidates.candidates[0].boundaries.clear();
        auto unspliced=project(mapped);
        require(tag(unspliced.records[0],"jM")=="c,-1" &&
            tag(unspliced.records[0],"jI")=="i,-1","unspliced tags differ");
        candidates.candidates[0].orientation=Orientation::Reverse;
        auto reverse=project(mapped);
        require(reverse.records[0].sequence=="TTTTTT" && reverse.records[0].quality=="FEDCBA","reverse output differs");
        config.outSAMmode="NoQS";
        require(project(mapped).records[0].quality=="*","NoQS differs");
        require(project({MappingClass::Unmapped,UnmappedReason::NoGoodWindow,0}).records[0].quality=="ABCDEF","unmapped quality changed");
        read.original.format=FastxFormat::Fasta;
        config.outSAMmode="Full";
        require(project(mapped).records[0].quality=="*","FASTA quality changed");
    } catch(const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
}
