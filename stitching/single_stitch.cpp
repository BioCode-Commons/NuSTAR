#include "single_stitch.hpp"
#include "stitch_detail.hpp"
#include "../reference/genome.hpp"

#include <algorithm>
#include <limits>

namespace nustar::star2711b {
namespace {

Failure stitchError(ErrorKind kind, std::string message) {
    return Failure{kind, Stage::Stitch, std::nullopt, 1, std::move(message)};
}

stitch_detail::Parameters parametersFor(const SingleStitchConfig&, const ScoreConfig&);

class StarStitchWorkspace final : public StageWorkspace {
public:
    StarStitchWorkspace(std::string reference_id, stitch_detail::Parameters parameters)
        : parameters(std::move(parameters)), identity_{"star2711b-single-stitch", 1, std::move(reference_id)} {}
    const WorkspaceIdentity& identity() const noexcept override { return identity_; }
    void reset() noexcept override {}
    void prepare(std::size_t seeds, std::size_t transcripts) {
        alignments.resize(seeds);
        if (included_capacity < seeds) {
            included = std::make_unique<bool[]>(seeds);
            included_capacity = seeds;
        }
        std::fill_n(included.get(), seeds, false);
        storage.resize(transcripts);
        ordered.resize(transcripts);
        for (std::size_t i = 0; i < transcripts; ++i)
            ordered[i] = &storage[i];
    }
    stitch_detail::Parameters parameters;
    std::vector<stitch_detail::uiWA> alignments;
    std::unique_ptr<bool[]> included;
    std::vector<stitch_detail::Transcript> storage;
    std::vector<stitch_detail::Transcript*> ordered;
private:
    WorkspaceIdentity identity_;
    std::size_t included_capacity = 0;
};

stitch_detail::Parameters parametersFor(const SingleStitchConfig& config,
                                         const ScoreConfig& score) {
    stitch_detail::Parameters result;
    const auto& common = config.common;
    result.alignIntronMin = common.alignIntronMin;
    result.alignIntronMax = common.alignIntronMax;
    result.alignSJoverhangMin = common.alignSJoverhangMin;
    result.alignSJDBoverhangMin = common.alignSJDBoverhangMin;
    result.alignSJstitchMismatchNmax = common.alignSJstitchMismatchNmax;
    result.alignSplicedMateMapLmin = common.alignSplicedMateMapLmin;
    result.alignSplicedMateMapLminOverLmate = common.alignSplicedMateMapLminOverLmate;
    result.alignTranscriptsPerWindowNmax = common.alignTranscriptsPerWindowNmax;
    result.alignSoftClipAtReferenceEnds.yes = common.alignSoftClipAtReferenceEnds;
    result.alignInsertionFlush.flushRight = common.insertion_flush_right;
    result.alignEndsType.ext[0][0] = config.force_end_extension[0];
    result.alignEndsType.ext[0][1] = config.force_end_extension[1];
    result.outFilterIntronMotifs = common.outFilterIntronMotifs;
    result.outFilterIntronStrands = common.outFilterIntronStrands;
    result.outSAMstrandField.type = common.outSAMstrandField == "intronMotif" ? 1 : 0;
    result.outFilterMismatchNoverLmax = common.outFilterMismatchNoverLmax;
    result.outFilterMultimapScoreRange = common.outFilterMultimapScoreRange;
    result.scoreGap = score.scoreGap;
    result.scoreGapNoncan = score.scoreGapNoncan;
    result.scoreGapGCAG = score.scoreGapGCAG;
    result.scoreGapATAC = score.scoreGapATAC;
    result.scoreDelOpen = score.scoreDelOpen;
    result.scoreDelBase = score.scoreDelBase;
    result.scoreInsOpen = score.scoreInsOpen;
    result.scoreInsBase = score.scoreInsBase;
    result.scoreStitchSJshift = score.scoreStitchSJshift;
    result.pGe.sjdbScore = score.sjdbScore;
    result.scoreGenomicLengthLog2scale = score.scoreGenomicLengthLog2scale;
    return result;
}

AlignmentCandidate toCandidate(const stitch_detail::Transcript& source,
    CandidateId id, const AlignmentWindow& window,
    const ReferenceMetadata& metadata) {
    AlignmentCandidate result{};
    result.id = id;
    result.window = window.id;
    result.chromosome = window.chromosome;
    result.orientation = window.orientation;
    result.original_read_orientation = source.roStr == 0 ?
        Orientation::Forward : Orientation::Reverse;
    result.single_mate = MateId::First;
    result.read_start = ReadOffset{source.rStart};
    result.original_read_start = ReadOffset{source.roStart};
    result.read_length = source.rLength;
    result.genome_start = GenomeOffset{source.gStart};
    result.genome_length = source.gLength;
    result.chromosome_start = metadata.chromosomes[window.chromosome.value].start.value;
    result.blocks.reserve(source.nExons);
    if (source.nExons) result.boundaries.reserve(source.nExons-1);
    for (unsigned i = 0; i < source.nExons; ++i) {
        const auto* exon = source.exons[i];
        const auto sjdb = exon[stitch_detail::EX_sjA] == std::numeric_limits<stitch_detail::uint>::max() ?
            std::optional<JunctionId>{} : std::optional<JunctionId>{JunctionId{exon[stitch_detail::EX_sjA]}};
        result.blocks.push_back(AlignmentBlock{ReadOffset{exon[stitch_detail::EX_R]},
            GenomeOffset{exon[stitch_detail::EX_G]}, exon[stitch_detail::EX_L], MateId::First, sjdb});
        if (i+1 < source.nExons) {
            const auto code = source.canonSJ[i];
            const auto kind = code == -2 ? BoundaryKind::Insertion :
                              code == -1 ? BoundaryKind::Deletion :
                              code == -3 ? BoundaryKind::MateGap : BoundaryKind::Intron;
            const auto motif = code >= 0 && code <= 6 ?
                static_cast<Motif>(code) : Motif::Noncanonical;
            result.boundaries.push_back(AlignmentBoundary{kind, motif,
                source.sjAnnot[i], static_cast<Strand>(source.sjStr[i]),
                source.shiftSJ[i][0], source.shiftSJ[i][1]});
        }
    }
    result.metrics = AlignmentMetrics{source.maxScore, source.nMatch,
        source.nMM, source.mappedLength, source.extendL, source.nGap,
        source.lGap, source.nDel, source.lDel, source.nIns, source.lIns,
        source.nUnique, source.nAnchor,
        {source.intronMotifs[0], source.intronMotifs[1], source.intronMotifs[2]},
        static_cast<Strand>(source.sjMotifStrand), source.sjYes};
    return result;
}

class StarSingleStitcher final : public ISingleStitcher {
public:
    explicit StarSingleStitcher(std::shared_ptr<const GenomeReferenceData> reference)
        : reference_(std::move(reference)), semantic_id_(reference_->semanticId()) {
        const auto genome = reference_->starGenome();
        static_assert(sizeof(BaseCode) == sizeof(char));
        genome_.G = reinterpret_cast<const char*>(genome.data());
        genome_.nGenome = genome.size();
        for (const auto& chromosome : reference_->storage().metadata.chromosomes) {
            genome_.chrStart.push_back(chromosome.start.value);
            genome_.chrLength.push_back(chromosome.length);
        }
        for (const auto& junction : reference_->junctions()) {
            genome_.sjdbStart.push_back(junction.lookup_intron.begin.value);
            genome_.sjdbEnd.push_back(junction.lookup_intron.end.value-1);
            genome_.sjdbMotif.push_back(static_cast<unsigned char>(junction.motif));
            genome_.sjdbShiftLeft.push_back(junction.shift_left);
            genome_.sjdbShiftRight.push_back(junction.shift_right);
            genome_.sjdbStrand.push_back(static_cast<unsigned char>(junction.strand));
        }
        genome_.sjdbN = reference_->junctions().size();
    }

    std::string_view backendId() const noexcept override {
        return "star2711b-single-stitch";
    }

    Result<std::unique_ptr<StageWorkspace>> createWorkspace(
        const SingleStitchConfig& config, const ScoreConfig& score,
        const ReferenceView& reference, const CompatibilityLimits& limits,
        const WorkspaceBudget&) const override {
        if (reference.semantic_id != semantic_id_ ||
            limits.max_blocks < stitch_detail::MAX_N_EXONS ||
            config.common.alignTranscriptsPerWindowNmax == 0)
            return stitchError(ErrorKind::Config, "single stitch reference or capacity mismatch");
        return std::unique_ptr<StageWorkspace>(
            std::make_unique<StarStitchWorkspace>(semantic_id_, parametersFor(config, score)));
    }

    Status stitchSingle(const SingleMappingRead& input,
        const WindowBatch& windows, const ReferenceView& reference,
        const SingleStitchConfig& config, const ScoreConfig&,
        StageWorkspace& workspace, CandidateSet& out) const override {
        const auto& read = input.view();
        out = CandidateSet{};
        out.generation = reference.generation;
        out.markers = windows.markers;
        auto& scratch = static_cast<StarStitchWorkspace&>(workspace);
        auto& parameters = scratch.parameters;
        stitch_detail::ReadAlign aligner;
        aligner.outFilterMismatchNmaxTotal = read.budget.mismatch_limit;
        aligner.readLength[0] = read.mates[0].clipped_length;
        std::uint64_t retained_transcripts = 0;
        for (const auto& window : windows.windows) {
            if (!window.alive || window.seeds.empty()) continue;
            if (parameters.alignTranscriptsPerWindowNmax >=
                    config.common.alignTranscriptsPerReadNmax ||
                retained_transcripts >= config.common.alignTranscriptsPerReadNmax -
                    parameters.alignTranscriptsPerWindowNmax) {
                out.markers.diagnostics.push_back(Diagnostic{
                    Diagnostic::Level::Warning, Stage::Stitch, {}, std::nullopt,
                    "alignTranscriptsPerReadNmax", "remaining stitch windows were skipped"});
                break;
            }
            scratch.prepare(window.seeds.size(),
                parameters.alignTranscriptsPerWindowNmax+1);
            auto& alignments = scratch.alignments;
            for (std::size_t i = 0; i < window.seeds.size(); ++i) {
                const auto& seed = window.seeds[i];
                alignments[i] = {seed.length, seed.read_start.value,
                    seed.genome_start.value, seed.repeat_count,
                    seed.anchor ? 1ULL : 0ULL,
                    static_cast<stitch_detail::uint>(seed.mate),
                    seed.sjdb ? seed.sjdb->value :
                        std::numeric_limits<stitch_detail::uint>::max()};
            }
            if (window.last_anchor && *window.last_anchor < alignments.size())
                alignments[*window.last_anchor][stitch_detail::WA_Anchor] = 2;
            stitch_detail::Transcript initial;
            initial.Chr = window.chromosome.value;
            initial.Str = static_cast<unsigned>(window.orientation);
            initial.roStr = read.revert_strand ? 1-initial.Str : initial.Str;
            const auto& oriented = read.encoded[initial.roStr == 0 ? 0 : 2];
            const auto* bases = reinterpret_cast<const char*>(oriented.data());
            auto& ordered = scratch.ordered;
            *ordered[0] = initial;
            stitch_detail::uint count = 0;
            stitch_detail::stitchWindowAligns(0, alignments.size(), 0,
                scratch.included.get(), 0, 0, initial, read.budget.legacy_Lread,
                alignments.data(), bases, genome_, parameters,
                ordered.data(), &count, &aligner);
            if (count == 0) continue;
            retained_transcripts += count;
            WindowCandidates window_candidates{window.id, {}};
            for (std::size_t i = 0; i < count; ++i) {
                const auto id = CandidateId{out.candidates.size()};
                out.candidates.push_back(toCandidate(*ordered[i], id, window,
                                                      *reference.metadata));
                window_candidates.ordered.push_back(id);
            }
            const auto first = window_candidates.ordered.front();
            if (!out.best ||
                out.candidates[first.value].metrics.score >
                    out.candidates[out.best->value].metrics.score ||
                (out.candidates[first.value].metrics.score ==
                    out.candidates[out.best->value].metrics.score &&
                 out.candidates[first.value].genome_length <
                    out.candidates[out.best->value].genome_length))
                out.best = first;
            out.windows.push_back(std::move(window_candidates));
        }
        if (!out.best || out.candidates[out.best->value].metrics.score == 0) {
            out.best.reset();
            out.windows.clear();
            out.candidates.clear();
            out.markers.legacy_map_marker = 999903;
        }
        return std::monostate{};
    }

private:
    std::shared_ptr<const GenomeReferenceData> reference_;
    std::string semantic_id_;
    stitch_detail::Genome genome_;
};

} // namespace

Result<std::shared_ptr<const ISingleStitcher>> makeStarSingleStitcher(
    ReferenceLease reference) {
    auto concrete = std::dynamic_pointer_cast<const GenomeReferenceData>(reference);
    if (!concrete) return stitchError(ErrorKind::Config, "STAR stitcher needs a NuSTAR reference");
    return std::shared_ptr<const ISingleStitcher>(
        std::make_shared<StarSingleStitcher>(std::move(concrete)));
}

} // namespace nustar::star2711b
