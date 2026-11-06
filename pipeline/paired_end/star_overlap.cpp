#include "overlap.hpp"
#include "../../reference/genome.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

// STAR 2.7.11b: ReadAlign_peOverlapMergeMap.cpp, SequenceFuns.cpp,
// and Transcript_alignScore.cpp, adapted to candidate IDs and module workspaces.
namespace nustar::star2711b {
namespace {

// STAR SequenceFuns.cpp::localSearchNisMM: N is a mismatch and the
// denominator is the number of matches, not the length of the overlap.
std::size_t overlapStart(std::span<const BaseCode> x, std::span<const BaseCode> y,
                         double mismatch_fraction) {
    std::size_t best=x.size(), best_matches=0, best_mismatches=0;
    for (std::size_t start=0;start<x.size();++start) {
        std::size_t matches=0, mismatches=0;
        for (std::size_t i=0;i<std::min(y.size(),x.size()-start);++i) {
            if (x[start+i]==y[i] && y[i]<=BaseCode::T) ++matches;
            else ++mismatches;
        }
        if ((matches>best_matches || (matches==best_matches && mismatches<best_mismatches)) &&
            double(mismatches)/double(matches)<=mismatch_fraction) {
            best=start; best_matches=matches; best_mismatches=mismatches;
        }
    }
    return best;
}

void rescore(AlignmentCandidate& candidate, const MappingRead& read,
             std::span<const BaseCode> genome, const ScoreConfig& config) {
    auto& metrics=candidate.metrics;
    metrics.score=0; metrics.matches=0; metrics.mismatches=0;
    const auto& bases=read.encoded[candidate.original_read_orientation==Orientation::Forward ? 0 : 2];
    for (const auto& block:candidate.blocks)
        for (std::size_t i=0;i<block.length;++i) {
            const auto r=bases[block.read_start.value+i], g=genome[block.genome_start.value+i];
            if (r>BaseCode::T || g>BaseCode::T) continue;
            if (r==g) { ++metrics.score; ++metrics.matches; }
            else { --metrics.score; ++metrics.mismatches; }
        }
    for (std::size_t i=0;i<candidate.boundaries.size();++i) {
        const auto& boundary=candidate.boundaries[i];
        const auto& a=candidate.blocks[i]; const auto& b=candidate.blocks[i+1];
        if (boundary.annotation==1) { metrics.score+=config.sjdbScore; continue; }
        switch (boundary.kind) {
        case BoundaryKind::MateGap: break;
        case BoundaryKind::Insertion:
            metrics.score+=Score(b.read_start.value-a.read_start.value-a.length)*config.scoreInsBase+config.scoreInsOpen;
            break;
        case BoundaryKind::Deletion:
            metrics.score+=Score(b.genome_start.value-a.genome_start.value-a.length)*config.scoreDelBase+config.scoreDelOpen;
            break;
        case BoundaryKind::Intron:
            metrics.score+=config.scoreGap;
            switch (boundary.motif) {
            case Motif::Noncanonical: metrics.score+=config.scoreGapNoncan; break;
            case Motif::GC_AG: case Motif::CT_GC: metrics.score+=config.scoreGapGCAG; break;
            case Motif::AT_AC: case Motif::GT_AT: metrics.score+=config.scoreGapATAC; break;
            default: break;
            }
            break;
        }
    }
    if (config.scoreGenomicLengthLog2scale!=0) {
        const auto& last=candidate.blocks.back();
        const auto length=std::max<std::uint64_t>(1,last.genome_start.value+last.length-candidate.blocks.front().genome_start.value);
        metrics.score+=Score(std::ceil(std::log2(double(length))*config.scoreGenomicLengthLog2scale-0.5));
    }
}

// Transcript::peOverlapSEtoPE, including its scalar-field assignments.
AlignmentCandidate project(const AlignmentCandidate& source, const MappingRead& pair,
    const std::array<std::uint64_t,2>& starts, std::uint64_t merged_length,
    std::uint64_t max_blocks) {
    AlignmentCandidate result{};
    result.id=source.id; result.window=source.window; result.chromosome=source.chromosome;
    result.orientation=source.orientation; result.original_read_orientation=source.original_read_orientation;
    result.genome_start=source.genome_start; result.genome_length=source.genome_length;
    result.chromosome_start=source.chromosome_start;
    result.metrics=source.metrics;
    result.metrics.extension_length=0;
    result.metrics.deletion_bases=source.metrics.deletion_count; // STAR assigns lDel=t.nDel here.
    result.metrics.has_splice=false;
    auto begin=starts;
    const auto strand=static_cast<unsigned>(source.orientation);
    const std::array<std::uint64_t,2> lengths{pair.mates[strand].clipped_length,pair.mates[1-strand].clipped_length};
    if (strand) {
        for (unsigned m=0;m<2;++m) begin[m]=merged_length-pair.mates[m].clipped_length-begin[m];
        std::swap(begin[0],begin[1]);
    }
    const AlignmentBoundary mate_gap{BoundaryKind::MateGap,Motif::Noncanonical,0,Strand::Unknown,0,0};
    for (unsigned mate=0;mate<2;++mate) {
        const auto end=begin[mate]+lengths[mate];
        for (std::size_t i=0;i<source.blocks.size();++i) {
            const auto& block=source.blocks[i];
            const auto left=std::max(begin[mate],block.read_start.value);
            const auto right=std::min(end,block.read_start.value+block.length);
            if (left>=right) continue;
            result.blocks.push_back({ReadOffset{left-begin[mate]+(mate ? lengths[0]+1 : 0)},
                GenomeOffset{block.genome_start.value+left-block.read_start.value},right-left,
                static_cast<MateId>(mate ? 1-strand : strand),block.sjdb});
            result.boundaries.push_back(i<source.boundaries.size() ? source.boundaries[i] : mate_gap);
            if (result.blocks.size()>max_blocks) { result.blocks.clear(); return result; }
        }
        if (!result.boundaries.empty()) result.boundaries.back()=mate_gap;
    }
    if (result.blocks.empty()) return result;
    result.boundaries.pop_back();
    for (const auto& block:result.blocks) result.read_length+=block.length;
    result.metrics.mapped_length=result.read_length;
    result.read_start=result.blocks.front().read_start;
    result.original_read_start=ReadOffset{result.original_read_orientation==Orientation::Forward ?
        result.read_start.value : pair.budget.legacy_Lread-result.read_start.value-result.read_length};
    return result;
}

class OverlapWorkspace final : public StageWorkspace {
public:
    explicit OverlapWorkspace(std::unique_ptr<StageWorkspace> mapper, std::uint64_t blocks)
        : mapper(std::move(mapper)), max_blocks(blocks) {}
    const WorkspaceIdentity& identity() const noexcept override { return identity_; }
    void reset() noexcept override { mapper->reset(); }
    std::unique_ptr<StageWorkspace> mapper;
    std::uint64_t max_blocks;
    CandidateSet merged_candidates;
private:
    WorkspaceIdentity identity_{"star2711b-paired-overlap",1,{}};
};

class StarPairedOverlap final : public IPairedOverlapRefiner {
public:
    explicit StarPairedOverlap(std::shared_ptr<const GenomeReferenceData> reference)
        : reference_(std::move(reference)) {}
    std::string_view backendId() const noexcept override { return "star2711b-paired-overlap"; }
    Result<std::unique_ptr<StageWorkspace>> createWorkspace(const OverlapConfig&,
        const ISingleCandidateMapper& mapper, const CompatibilityLimits& limits,
        const WorkspaceBudget& budget) const override {
        auto made=mapper.createWorkspace(budget);
        if (auto* failure=std::get_if<Failure>(&made)) return *failure;
        return std::unique_ptr<StageWorkspace>(std::make_unique<OverlapWorkspace>(
            std::get<std::unique_ptr<StageWorkspace>>(std::move(made)),limits.max_blocks));
    }
    Status refinePair(const PairedMappingRead& input, const ReferenceView&,
        const OverlapConfig& config, const ScoreConfig& score, const ISingleCandidateMapper& mapper,
        StageWorkspace& workspace, CandidateSet& candidates, OverlapState& out) const override {
        out={}; out.attempted=true;
        const auto& pair=input.view();
        const auto n0=pair.mates[0].clipped_length,n1=pair.mates[1].clipped_length;
        const auto bases=std::span<const BaseCode>(pair.encoded[0]);
        const auto first=bases.first(n0),second=bases.subspan(n0+1,n1);
        const auto s1=overlapStart(first,second,config.peOverlapMMp);
        const auto s0=overlapStart(second,first,config.peOverlapMMp);
        const auto o1=std::min(n1,n0-s1),o0=std::min(n0,n1-s0);
        const auto overlap=std::max(o0,o1);
        if (overlap<config.peOverlapNbasesMin || overlap==0) return std::monostate{};
        out.overlapping_bases=overlap;
        out.mate_starts_in_merged=o1>=o0 ? std::array<std::uint64_t,2>{0,s1} : std::array<std::uint64_t,2>{s0,0};
        MappingRead merged{};
        merged.fragment=pair.fragment; merged.revert_strand=pair.revert_strand;
        auto& sequence=merged.encoded[0];
        sequence.reserve(n0+n1-overlap);
        const auto leading=o1>=o0 ? first : second, trailing=o1>=o0 ? second : first;
        sequence.insert(sequence.end(),leading.begin(),leading.end());
        sequence.insert(sequence.end(),trailing.begin()+overlap,trailing.end());
        const auto length=sequence.size();
        merged.mates.push_back({length,length,{}, {ReadOffset{0},ReadOffset{length}},false});
        merged.budget={pair.budget.mismatch_limit,length,pair.budget.legacy_original_pair_length,length};
        for (auto b:sequence) merged.encoded[1].push_back(b<=BaseCode::T ? static_cast<BaseCode>(static_cast<unsigned>(b)^3) : b);
        merged.encoded[2].assign(merged.encoded[1].rbegin(),merged.encoded[1].rend());
        auto& scratch=static_cast<OverlapWorkspace&>(workspace);
        auto status=mapper.mapSingle(SingleMappingRead::fromPrepared(std::move(merged)),*scratch.mapper,scratch.merged_candidates);
        if (std::holds_alternative<Failure>(status)) return status;
        auto& mapped=scratch.merged_candidates;
        if (mapped.windows.empty()) return std::monostate{};
        const auto old_score=candidates.best ? candidates.candidates[candidates.best->value].metrics.score : Score(0);
        CandidateSet converted{};
        converted.generation=candidates.generation; converted.markers=candidates.markers;
        for (const auto& window:mapped.windows) {
            WindowCandidates converted_window{window.window,{}};
            for (const auto id:window.ordered) {
                auto candidate=project(mapped.candidates[id.value],pair,out.mate_starts_in_merged,length,scratch.max_blocks);
                if (candidate.blocks.empty()) continue;
                rescore(candidate,pair,reference_->starGenome(),score);
                candidate.id=CandidateId{converted.candidates.size()};
                converted_window.ordered.push_back(candidate.id);
                converted.candidates.push_back(std::move(candidate));
                auto& ordered=converted_window.ordered;
                if (converted.candidates[ordered.back().value].metrics.score>converted.candidates[ordered.front().value].metrics.score)
                    std::swap(ordered.back(),ordered.front());
            }
            if (converted_window.ordered.empty()) continue;
            const auto best=converted_window.ordered.front();
            if (!converted.best || converted.candidates[best.value].metrics.score>converted.candidates[converted.best->value].metrics.score)
                converted.best=best;
            converted.windows.push_back(std::move(converted_window));
        }
        // STAR replaces the candidates even when the new score is lower;
        // peOv.yes records the score comparison, not whether to restore PE.
        out.accepted=converted.best && old_score<=converted.candidates[converted.best->value].metrics.score;
        candidates=std::move(converted);
        return std::monostate{};
    }
private:
    std::shared_ptr<const GenomeReferenceData> reference_;
};
}
Result<std::shared_ptr<const IPairedOverlapRefiner>> makeStarPairedOverlapRefiner(ReferenceLease reference) {
    auto concrete=std::dynamic_pointer_cast<const GenomeReferenceData>(reference);
    if (!concrete) return Failure{ErrorKind::Config,Stage::Plan,{},1,"paired overlap needs a genome reference"};
    return std::shared_ptr<const IPairedOverlapRefiner>(std::make_shared<StarPairedOverlap>(std::move(concrete)));
}
} // namespace nustar::star2711b
