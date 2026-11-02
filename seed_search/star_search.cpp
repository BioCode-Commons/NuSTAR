#include "star_search.hpp"
#include "base_compare.hpp"

#include <algorithm>
#include <atomic>
#include <array>
#include <cmath>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <utility>

namespace nustar::star2711b {
namespace {

Failure seedError(ErrorKind kind, std::string message,
                  std::optional<FragmentId> fragment = {}) {
    return Failure{kind, Stage::Seed, fragment, 1, std::move(message)};
}

std::uint64_t next_source_id() {
    static std::atomic<std::uint64_t> next{1};
    // Opaque IDs need uniqueness, not a process-wide per-read ordering.
    // Amortize the shared cache-line update across each worker's local block.
    thread_local std::uint64_t current = 0, end = 0;
    if (current == end) {
        current = next.fetch_add(4096, std::memory_order_relaxed);
        end = current+4096;
    }
    return current++;
}

} // namespace

StarSeedSearch::StarSeedSearch(std::shared_ptr<const SuffixArrayIndex> index,
    std::shared_ptr<const GenomeReferenceData> reference, SeedConfig policy)
    : index_(std::move(index)), reference_(std::move(reference)),
      policy_(policy) {
    if (!index_ || !reference_)
        throw std::invalid_argument("seed search requires index and reference");
    real_length_ = reference_->storage().forward.size();
    if (index_->identity().reference_semantic_id != reference_->semanticId() ||
        index_->physicalLength() != 2*(real_length_+reference_->synthetic().size())+80)
        throw std::invalid_argument("seed search index/reference identity mismatch");
    if (policy_.seedSearchStartLmaxOverLread < 0 ||
        !std::isfinite(policy_.seedSearchStartLmaxOverLread) ||
        policy_.seedSplitMin == 0 || policy_.seedMapMin == 0 ||
        policy_.maxNsplit == 0 || policy_.seedPerReadNmax == 0)
        throw std::invalid_argument("invalid STAR seed policy");
    genome_ = reference_->starGenome();

}

// STAR SuffixArrayFuns.cpp::compareSeqToGenome. The four direction cases
// use the same shared Genome.G bytes. Only legacy physical SA decoding differs;
// suffix ranks and all comparison decisions are unchanged.
template<bool WordComparison>
std::uint64_t StarSeedSearch::compare(const MappingRead& read, std::uint64_t start,
    std::uint64_t length, std::uint64_t known, std::uint64_t rank,
    bool forward, bool& query_greater) const {
    if (known == length) return length;
    const auto position = index_->suffixArray()[rank];
    const auto synthetic_start = 2*real_length_+40;
    const auto reverse_synthetic_start = synthetic_start+reference_->synthetic().size()+20;
    const BaseCode* g;
    bool genome_forward;
    if (position < real_length_) {
        g = genome_.data()+position;
        genome_forward = true;
    } else if (position < synthetic_start) {
        g = genome_.data()+(2*real_length_+19-position);
        genome_forward = false;
    } else if (position < reverse_synthetic_start) {
        g = genome_.data()+real_length_+position-synthetic_start;
        genome_forward = true;
    } else {
        g = genome_.data()+genome_.size()-1-(position-reverse_synthetic_start);
        genome_forward = false;
    }
    const auto* s = read.encoded[forward == genome_forward ? 0 : 1].data()+start;
    if constexpr (WordComparison) {
        s += forward ? static_cast<std::ptrdiff_t>(known) : -static_cast<std::ptrdiff_t>(known);
        g += genome_forward ? static_cast<std::ptrdiff_t>(known) : -static_cast<std::ptrdiff_t>(known);
        const auto remaining=length-known;
        const auto* begin=genome_.data();
        const auto* end=begin+genome_.size();
        if (forward && genome_forward)
            return known+detail::compareBases<true,true>(s,g,remaining,begin,end,query_greater);
        if (forward)
            return known+detail::compareBases<true,false>(s,g,remaining,begin,end,query_greater);
        if (genome_forward)
            return known+detail::compareBases<false,true>(s,g,remaining,begin,end,query_greater);
        return known+detail::compareBases<false,false>(s,g,remaining,begin,end,query_greater);
    } else {
        const auto remaining = length-known;
        if (forward && genome_forward) {
            s += known; g += known;
            for (std::uint64_t i=0; i<remaining; ++i)
                if (s[i] != g[i]) { query_greater = s[i] > g[i]; return i+known; }
        } else if (forward && !genome_forward) {
            s += known; g -= known;
            for (std::int64_t i=0; std::uint64_t(i)<remaining; ++i)
                if (s[i] != g[-i]) {
                    query_greater = !(s[i] > g[-i] || std::uint8_t(g[-i]) > 3);
                    return i+known;
                }
        } else if (!forward && genome_forward) {
            s -= known; g += known;
            for (std::int64_t i=0; std::uint64_t(i)<remaining; ++i)
                if (s[-i] != g[i]) { query_greater = s[-i] > g[i]; return i+known; }
        } else {
            s -= known; g -= known;
            for (std::int64_t i=0; std::uint64_t(i)<remaining; ++i)
                if (s[-i] != g[-i]) {
                    query_greater = !(s[-i] > g[-i] || std::uint8_t(g[-i]) > 3);
                    return i+known;
                }
        }
        return length;
    }
}

StarSeedSearch::Match StarSeedSearch::longest(const MappingRead& read,
    std::uint64_t start, std::uint64_t length, bool forward) const {
    const auto& suffixes = index_->suffixArray();
    if (suffixes.empty() || length == 0) return {};
    const auto query = forward ? std::span(read.encoded[0]).subspan(start, length) :
        std::span(read.encoded[2]).subspan(read.budget.legacy_Lread-1-start, length);
    std::size_t domain_begin = 0, domain_end = suffixes.size();
    std::uint32_t known = 0;
    const auto& sai=index_->prefixIndex();
    unsigned level=std::min<std::uint64_t>(length,index_->prefixBases());
    std::uint64_t code=0;
    for (unsigned i=0;i<level;++i) code=(code<<2)|std::uint8_t(query[i]);
    std::uint64_t first=0;
    while (level) {
        first=sai[sai.levelStart(level)+code];
        if (!(first&sai.absentMask())) break;
        --level; code>>=2;
    }
    // STAR assumes a present one-base prefix. Empty/degenerate references
    // retain the existing full-domain path instead of indexing level -1.
    if (level) {
        bool upper_good=false;
        std::uint64_t last=suffixes.size()-1;
        if (sai.levelStart(level)+code+1<sai.levelEnd(level)) {
            const auto next=sai[sai.levelStart(level)+code+1];
            if (!(next&sai.absentMask())) {
                last=(next&~sai.nMask())-1;
                upper_good=true;
            }
        }
        const bool no_n=!(first&sai.nMask());
        domain_begin=first&sai.rankMask(); domain_end=last+1;
        if (level<index_->prefixBases() && no_n && upper_good)
            return {level,domain_begin,last,domain_end-domain_begin};
        if (first==last && no_n && upper_good) {
            bool greater=false;
            // Only the existing singleton-extension branch uses word loads.
            // Interval probes usually differ early and are faster bytewise.
            const auto matched=compare<true>(read,start,length,level,first,forward,greater);
            return {matched,first,first,1};
        }
        if (upper_good && no_n) known=level;
    }
    // STAR SuffixArrayFuns.cpp::maxMappableLength and findMultRange.
    // Keep endpoint match lengths and their two preceding positions, including
    // STAR's exact-match exit and its tie choice, including SAi fallback bounds.
    bool greater = false;
    auto cmp = [&](std::uint64_t rank, std::uint64_t n, std::uint64_t equal) {
        return compare<false>(read, start, n, equal, rank, forward, greater);
    };
    std::uint64_t i1=domain_begin, i2=domain_end-1;
    auto L1=cmp(i1,length,known), L2=cmp(i2,length,known);
    auto L=std::min(L1,L2);
    auto L1a=L1, L1b=L1, i1a=i1, i1b=i1;
    auto L2a=L2, L2b=L2, i2a=i2, i2b=i2;
    auto i3=i1, L3=L1;
    while (i1+1<i2) {
        i3=i1+(i2-i1)/2;
        L3=cmp(i3,length,L);
        if (L3==length) break;
        if (greater) {
            if (L3>L1) { L1b=L1a; L1a=L1; i1b=i1a; i1a=i1; }
            i1=i3; L1=L3;
        } else {
            if (L3>L2) { L2b=L2a; L2a=L2; i2b=i2a; i2a=i2; }
            i2=i3; L2=L3;
        }
        L=std::min(L1,L2);
    }
    if (L3<length) {
        if (L1>L2) { i3=i1; L3=L1; }
        else { i3=i2; L3=L2; }
    }
    auto boundary = [&](std::uint64_t i, std::uint64_t len,
                        std::uint64_t ia, std::uint64_t la,
                        std::uint64_t ib, std::uint64_t lb) {
        if (len<L3) { lb=len; ib=i; ia=i3; }
        else if (la<len) { lb=la; ib=ia; ia=i; }
        while (ib+1<ia || ib>ia+1) {
            const auto ic=ia/2+ib/2+(ia%2+ib%2)/2;
            const auto lc=cmp(ic,L3,lb);
            if (lc==L3) ia=ic;
            else { ib=ic; lb=lc; }
        }
        return ia;
    };
    const auto begin=boundary(i1,L1,i1a,L1a,i1b,L1b);
    const auto end=boundary(i2,L2,i2a,L2a,i2b,L2b);
    return {L3,begin,end,end-begin+1};
}

Result<StarSeedSearchResult> StarSeedSearch::search(const MappingRead& read) const {
    StarSeedSearchResult result;
    if (auto status = search(read, result); std::holds_alternative<Failure>(status))
        return std::get<Failure>(std::move(status));
    return result;
}

Status StarSeedSearch::search(const MappingRead& read, StarSeedSearchResult& result) const {
    result.splits.clear();
    result.pieces.clear();
    result.counters = {};
    result.map_marker = 0;
    const auto length = read.budget.legacy_Lread;
    if (length == 0) return std::monostate{};
    auto& splits = result.splits;
    splits.reserve(2);
    result.pieces.reserve(32);
    std::uint64_t position = 0;
    std::uint64_t mate = 0;
    while (position < length && splits.size() < policy_.maxNsplit) {
        while (position < length && static_cast<std::uint8_t>(read.encoded[0][position]) > 3) {
            if (read.encoded[0][position] == BaseCode::MateSpacer) ++mate;
            ++position;
        }
        if (position == length) break;
        const auto begin = position;
        while (position < length &&
               static_cast<std::uint8_t>(read.encoded[0][position]) <= 3) ++position;
        const auto piece_length = position-begin;
        result.counters.no_good_piece_length =
            std::max(result.counters.no_good_piece_length, piece_length);
        if (piece_length >= policy_.seedSplitMin)
            splits.push_back(SplitSegment{
                ReadRange{ReadOffset{begin}, ReadOffset{position}},
                mate == 0 ? MateId::First : MateId::Second});
    }

    auto store = [&](std::uint64_t direction, std::uint64_t shift,
                     Match match, std::uint64_t fragment) -> Status {
        auto& count = result.counters;
        if (match.multiplicity > policy_.seedMultimapNmax) {
            if (count.min_multimaps == 0 || match.multiplicity < count.min_multimaps) {
                count.min_multimaps = match.multiplicity;
                count.min_multimap_length = match.length;
            }
            return std::monostate{};
        }
        count.all_alignments += match.multiplicity;
        (match.multiplicity == 1 ? count.unique_alignments :
                                   count.multiple_alignments) += match.multiplicity;
        const auto read_start = direction == 0 ? shift : shift+1-match.length;
        auto& pieces = result.pieces;
        // STAR storeAligns scans backward: equal starts keep longer pieces
        // first, and an equal (start,length) is discarded before insertion.
        auto place = pieces.end();
        while (place != pieces.begin()) {
            const auto& previous = *(place-1);
            if (previous.read_start <= read_start) {
                if (previous.read_start != read_start || previous.length >= match.length) {
                    if (previous.read_start == read_start && previous.length == match.length)
                        return std::monostate{};
                    break;
                }
            }
            --place;
        }
        if (pieces.size() >= policy_.seedPerReadNmax)
            return seedError(ErrorKind::Capacity, "STAR seedPerReadNmax exceeded",
                             read.fragment);
        pieces.insert(place, StarSeedPiece{read_start, match.length, direction,
            match.multiplicity, match.begin, match.end, fragment});
        if (match.multiplicity == 1)
            count.max_unique_length = std::max(count.max_unique_length, match.length);
        else {
            if (count.min_multimaps == 0 || match.multiplicity < count.min_multimaps) {
                count.min_multimaps = match.multiplicity;
                count.min_multimap_length = match.length;
            }
        }
        return std::monostate{};
    };

    const auto effective_start = std::min(policy_.seedSearchStartLmax,
        static_cast<std::uint64_t>(policy_.seedSearchStartLmaxOverLread*(length-1)));
    for (const auto& split : splits) {
        const auto begin = split.read.begin.value;
        const auto split_length = split.read.end.value-begin;
        const auto fragment = split.mate == MateId::First ? 0u : 1u;
        if (policy_.seedSearchStartLmax > 0 && effective_start == 0 &&
            split_length > effective_start)
            return seedError(ErrorKind::Config, "STAR seed start step is zero",
                             read.fragment);
        const auto starts = policy_.seedSearchStartLmax > 0 &&
            effective_start < split_length ? split_length/effective_start+1 : 1;
        const auto step = split_length/starts;
        bool reverse_first_needed = true;
        for (std::uint64_t direction = 0; direction < 2; ++direction)
            for (std::uint64_t start = 0; start < starts; ++start) {
                if (reverse_first_needed || start > 0) {
                    std::uint64_t mapped = 0;
                    while (start*step+mapped+policy_.seedMapMin < split_length) {
                        const auto shift = direction == 0 ?
                            begin+start*step+mapped :
                            begin+split_length-start*step-1-mapped;
                        const auto remaining = split_length-mapped-start*step;
                        const auto match = longest(read, shift, remaining, direction == 0);
                        if (match.length == 0)
                            return seedError(ErrorKind::Index,
                                "zero-length STAR seed would not advance", read.fragment);
                        auto status = store(direction, shift, match, fragment);
                        if (std::holds_alternative<Failure>(status))
                            return std::get<Failure>(std::move(status));
                        // Preserve STAR's comparison with split length, including its
                        // offset behavior for segments after an N or mate spacer.
                        if (direction == 0 && start == 0 && mapped == 0 &&
                            shift+match.length == split_length)
                            reverse_first_needed = false;
                        mapped += match.length;
                    }
                }
                if (policy_.seedSearchLmax > 0) {
                    const auto shift = direction == 0 ? begin+start*step :
                        begin+split_length-start*step-1;
                    const auto available = direction == 0 ?
                        begin+split_length-shift : shift+1;
                    const auto requested = std::min(policy_.seedSearchLmax, available);
                    auto status = store(direction, shift,
                        longest(read, shift, requested, direction == 0), fragment);
                    if (std::holds_alternative<Failure>(status))
                        return std::get<Failure>(std::move(status));
                }
            }
    }
    if (length < policy_.outFilterMatchNmin) result.map_marker = 999910;
    else if (splits.empty()) result.map_marker = 999904;
    else if (result.counters.all_alignments == 0) result.map_marker = 999901;
    return std::monostate{};
}

namespace {

struct StarSeedResources {
    std::shared_ptr<const SuffixArrayIndex> index;
    std::shared_ptr<const GenomeReferenceData> reference;
    std::shared_ptr<const std::vector<std::uint32_t>> chromosome_bins;
};

class StarSeedOccurrenceSource;

class StarSeedHitCursor final : public ISeedHitCursor {
public:
    StarSeedHitCursor(std::shared_ptr<const StarSeedOccurrenceSource> source,
                      std::size_t piece_index);
    Status next(std::uint64_t max_hits, SeedHitPage& out) override;
private:
    std::shared_ptr<const StarSeedOccurrenceSource> source_;
    std::size_t piece_index_;
    std::uint64_t rank_;
    std::uint64_t emitted_ = 0;
};

class StarSeedOccurrenceSource final : public ISeedOccurrenceSource,
                                       public std::enable_shared_from_this<StarSeedOccurrenceSource> {
public:
    struct ProjectedHit {
        ChromosomeId chromosome;
        Orientation orientation;
        std::array<SeedHitBlock, 2> blocks;
        std::uint8_t block_count;
    };
    StarSeedOccurrenceSource(std::shared_ptr<const StarSeedResources> resources,
        SeedDeliveryContract delivery)
        : resources_(std::move(resources)), index_(resources_->index.get()),
          reference_(resources_->reference.get()),
          chromosome_bins_(resources_->chromosome_bins.get()) {
        identity_ = OccurrenceSourceIdentity{
            OccurrenceSourceId{}, {},
            reference_->view().generation, std::string(reference_->semanticId()),
            std::move(delivery)};
    }
    // Called only before publication. Published sources and cursors remain
    // immutable, including after the workspace is destroyed.
    void reset(const MappingRead& read) {
        identity_.source = OccurrenceSourceId{next_source_id()};
        identity_.fragment = read.fragment;
        read_length_ = read.budget.legacy_Lread;
        revert_strand_ = read.revert_strand;
    }
    StarSeedSearchResult native;
    const OccurrenceSourceIdentity& identity() const noexcept override { return identity_; }
    Result<std::unique_ptr<ISeedHitCursor>> open(OccurrenceHandle handle) const override {
        return std::unique_ptr<ISeedHitCursor>(
            std::make_unique<StarSeedHitCursor>(shared_from_this(), handle.set.value));
    }
    Status visit(OccurrenceHandle handle,
        const std::function<bool(const SeedHitView&)>& callback) const override {
        const auto& selected = native.pieces[handle.set.value];
        std::uint64_t emitted = 0;
        for (auto rank = selected.sa_begin; rank <= selected.sa_end; ++rank) {
            ProjectedHit hit;
            project(handle.set.value, rank, hit);
            if (hit.block_count && !callback(SeedHitView{
                    hit.chromosome, hit.orientation, emitted++,
                    std::span(hit.blocks).first(hit.block_count)}))
                break;
        }
        return std::monostate{};
    }
    Status visitBlocks(OccurrenceHandle handle,
        const std::function<bool(const SeedHitBlock&, Orientation)>& callback) const override {
        // STAR copies PC fields before walking the SA interval. Value locals
        // stay valid across the callback without reloading source metadata.
        const auto piece = native.pieces[handle.set.value];
        const auto real = reference_->storage().forward.size();
        const auto read_length = read_length_;
        const auto revert = revert_strand_;
        for (auto rank = piece.sa_begin; rank <= piece.sa_end; ++rank) {
            const auto [start, strand] = locate(piece, rank);
            const auto orientation = Orientation(strand ^ piece.direction ^ revert);
            const auto read_start = piece.direction != strand ?
                read_length-(piece.length+piece.read_start) : piece.read_start;
            if (start < real) {
                const SeedHitBlock block{ReadOffset{read_start}, GenomeOffset{start},
                    piece.length, piece.fragment == 0 ? MateId::First : MateId::Second, {}};
                if (!callback(block, orientation)) break;
            } else {
                ProjectedHit hit;
                hit.block_count = 0;
                projectJunction(piece, start, read_start, hit);
                for (unsigned i = 0; i < hit.block_count; ++i)
                    if (!callback(hit.blocks[i], orientation)) return std::monostate{};
            }
        }
        return std::monostate{};
    }
    Status visitAnchors(OccurrenceHandle handle,
        const std::function<bool(GenomeOffset, Orientation)>& callback) const override {
        const auto piece = native.pieces[handle.set.value];
        const auto real = reference_->storage().forward.size();
        for (auto rank = piece.sa_begin; rank <= piece.sa_end; ++rank) {
            const auto [start, strand] = locate(piece, rank);
            const auto orientation = Orientation(strand ^ piece.direction ^ revert_strand_);
            if (start < real) {
                if (!callback(GenomeOffset{start}, orientation)) break;
            } else {
                ProjectedHit hit;
                hit.block_count = 0;
                projectJunction(piece, start, 0, hit);
                for (unsigned i = 0; i < hit.block_count; ++i)
                    if (!callback(hit.blocks[i].genome_start, orientation))
                        return std::monostate{};
            }
        }
        return std::monostate{};
    }
    Status visitBatch(OccurrenceHandle handle,
        const std::function<bool(std::span<const SeedHitView>)>& callback) const override {
        const auto& selected=native.pieces[handle.set.value];
        std::array<ProjectedHit,16> storage;
        std::array<SeedHitView,16> views;
        std::size_t count=0;
        std::uint64_t emitted=0;
        for (auto rank=selected.sa_begin;rank<=selected.sa_end;++rank) {
            auto& hit=storage[count];
            project(handle.set.value,rank,hit);
            if (!hit.block_count) continue;
            views[count++]={hit.chromosome,hit.orientation,emitted++,
                std::span(hit.blocks).first(hit.block_count)};
            if (count==views.size()) {
                if (!callback(views)) return std::monostate{};
                count=0;
            }
        }
        if (count) callback(std::span(views).first(count));
        return std::monostate{};
    }
    const StarSeedPiece& piece(std::size_t index) const { return native.pieces[index]; }
    std::span<const StarSeedPiece> pieces() const { return native.pieces; }
    std::uint64_t suffix(std::uint64_t rank) const { return index_->suffixArray()[rank]; }
    // As in STAR's stitchPieces/sjAlignSplit, valid SA hits and chromosome
    // tables are supplied by the bound index; no per-occurrence range checks.
    [[gnu::always_inline]] inline void appendBlock(const StarSeedPiece& piece,
        ProjectedHit& hit, std::uint64_t genome_start, std::uint64_t length,
        std::uint64_t read_offset, std::optional<JunctionId> junction) const {
        const auto& storage=reference_->storage();
        const auto bin=genome_start >> storage.chromosome_bin_bits;
        hit.chromosome=storage.metadata.chromosomes[(*chromosome_bins_)[bin]].id;
        hit.blocks[hit.block_count++]=SeedHitBlock{ReadOffset{read_offset},
            GenomeOffset{genome_start},length,
            piece.fragment==0 ? MateId::First : MateId::Second,junction};
    }
    [[gnu::noinline]] void projectJunction(const StarSeedPiece& piece,
        std::uint64_t forward_start, std::uint64_t read_start, ProjectedHit& hit) const {
        const auto overhang=reference_->sjdbOverhang();
        const auto synthetic_offset=forward_start-reference_->storage().forward.size();
        const auto stride=2*overhang+1;
        const auto within=synthetic_offset%stride;
        // This is STAR's alignment decision, not defensive validation.
        if (within>=overhang || within+piece.length<=overhang) return;
        const auto& junction=reference_->junctions()[synthetic_offset/stride];
        const auto donor_length=overhang-within;
        appendBlock(piece,hit,junction.donor_sequence_begin.value+within,
                    donor_length,read_start,junction.id);
        appendBlock(piece,hit,junction.acceptor_sequence_begin.value,
                    piece.length-donor_length,read_start+donor_length,junction.id);
    }
    [[gnu::always_inline]] inline std::pair<std::uint64_t, std::uint8_t>
    locate(const StarSeedPiece& piece, std::uint64_t rank) const {
        const auto real = reference_->storage().forward.size();
        const auto synthetic = reference_->synthetic().size();
        const auto combined = real+synthetic;
        const auto position = suffix(rank);
        const auto forward_synthetic = 2*real+40;
        const auto reverse_synthetic = forward_synthetic+synthetic+20;
        std::uint64_t sa_offset = 0;
        std::uint8_t sa_strand = 0;
        if (position < real) sa_offset = position;
        else if (position < 2*real+20) {
            sa_offset = synthetic+position-(real+20);
            sa_strand = 1;
        } else if (position < reverse_synthetic) {
            sa_offset = real+position-forward_synthetic;
        } else {
            sa_offset = position-reverse_synthetic;
            sa_strand = 1;
        }
        const auto forward_start = sa_strand == 0 ? sa_offset :
            combined-(piece.length+sa_offset);
        return {forward_start, sa_strand};
    }
    [[gnu::always_inline]] inline void project(std::size_t piece_index, std::uint64_t rank, ProjectedHit& hit) const {
        const auto& piece = native.pieces[piece_index];
        const auto [forward_start, sa_strand] = locate(piece, rank);
        const auto real = reference_->storage().forward.size();
        const auto reverse_read_coordinates = piece.direction != sa_strand;
        const auto read_start = reverse_read_coordinates ?
            read_length_-(piece.length+piece.read_start) : piece.read_start;
        auto orientation = sa_strand ^ piece.direction;
        if (revert_strand_) orientation ^= 1;
        hit.chromosome = ChromosomeId{0};
        hit.orientation = orientation == 0 ? Orientation::Forward : Orientation::Reverse;
        hit.block_count = 0;
        if (forward_start>=real) return projectJunction(piece,forward_start,read_start,hit);
        appendBlock(piece,hit,forward_start,piece.length,read_start,{});
    }
private:
    friend class StarSeedHitCursor;
    // Each workspace owns a separate lease control block. Sources/cursors can
    // outlive that workspace, without every worker contending on the global
    // index/reference shared_ptr reference counters for every read.
    std::shared_ptr<const StarSeedResources> resources_;
    const SuffixArrayIndex* index_;
    const GenomeReferenceData* reference_;
    const std::vector<std::uint32_t>* chromosome_bins_;
    std::uint64_t read_length_;
    bool revert_strand_;
    OccurrenceSourceIdentity identity_;
};

StarSeedHitCursor::StarSeedHitCursor(
    std::shared_ptr<const StarSeedOccurrenceSource> source,
    std::size_t piece_index)
    : source_(std::move(source)), piece_index_(piece_index),
      rank_(source_->piece(piece_index).sa_begin) {}

Status StarSeedHitCursor::next(std::uint64_t max_hits, SeedHitPage& out) {
    if (max_hits == 0)
        return seedError(ErrorKind::Config, "seed hit page size must be positive");
    out.hits.clear();
    const auto& piece = source_->piece(piece_index_);
    while (rank_ <= piece.sa_end && out.hits.size() < max_hits) {
        StarSeedOccurrenceSource::ProjectedHit projected;
        source_->project(piece_index_, rank_++, projected);
        if (projected.block_count) {
            SeedHit owned{projected.chromosome, projected.orientation,
                          emitted_++, {}};
            owned.blocks.assign(projected.blocks.begin(),
                                projected.blocks.begin()+projected.block_count);
            out.hits.push_back(std::move(owned));
        }
    }
    out.end = rank_ > piece.sa_end;
    return std::monostate{};
}

// Recycle only after the last strong lease has been released. Each publication
// gets a fresh shared_ptr control block: old weak_ptrs stay expired, and cursors
// retained by another thread cannot observe the next read's PC/identity.
class StarSeedSourcePool : public std::enable_shared_from_this<StarSeedSourcePool> {
public:
    std::shared_ptr<StarSeedOccurrenceSource> acquire(
        const std::shared_ptr<const StarSeedResources>& resources,
        const SeedDeliveryContract& delivery) {
        std::unique_ptr<StarSeedOccurrenceSource> source;
        {
            std::lock_guard lock(mutex_);
            source = std::move(cached_);
        }
        if (!source) source = std::make_unique<StarSeedOccurrenceSource>(resources, delivery);
        return {source.release(), [pool = weak_from_this()](StarSeedOccurrenceSource* value) {
            std::unique_ptr<StarSeedOccurrenceSource> released(value);
            if (auto owner = pool.lock()) {
                std::lock_guard lock(owner->mutex_);
                if (!owner->cached_) owner->cached_ = std::move(released);
            }
        }};
    }
private:
    std::mutex mutex_;
    std::unique_ptr<StarSeedOccurrenceSource> cached_;
};

class StarSeedWorkspace final : public StageWorkspace {
public:
    StarSeedWorkspace(WorkspaceIdentity identity, StarSeedResources resources)
        : resources(std::make_shared<const StarSeedResources>(std::move(resources))),
          identity_(std::move(identity)) {}
    const std::shared_ptr<const StarSeedResources> resources;
    const std::shared_ptr<StarSeedSourcePool> sources = std::make_shared<StarSeedSourcePool>();
    const WorkspaceIdentity& identity() const noexcept override { return identity_; }
    void reset() noexcept override {}
private:
    WorkspaceIdentity identity_;
};

class StarSeedFinder final : public ISeedFinder {
public:
    StarSeedFinder(std::shared_ptr<const SuffixArrayIndex> index,
        std::shared_ptr<const GenomeReferenceData> reference,
        std::shared_ptr<const SeedSearchPlan> plan)
        : index_(std::move(index)), reference_(std::move(reference)),
          plan_(std::move(plan)),
          search_(index_, reference_, plan_->compatibility_policy),
          workspace_identity_{"star-2.7.11b-seed-v1", 1,
              index_->identity().artifact_sha256+":"+
              std::string(reference_->semanticId())} {
        const auto& storage = reference_->storage();
        const auto& chromosomes = storage.metadata.chromosomes;
        auto bins = std::make_shared<std::vector<std::uint32_t>>(
            (storage.forward.size() >> storage.chromosome_bin_bits)+1);
        std::uint32_t chromosome = 0;
        for (std::size_t bin = 0; bin < bins->size(); ++bin) {
            const auto position = std::uint64_t(bin) << storage.chromosome_bin_bits;
            while (chromosome+1 < chromosomes.size() &&
                   position >= chromosomes[chromosome+1].start.value)
                ++chromosome;
            (*bins)[bin] = chromosome;
        }
        chromosome_bins_ = std::move(bins);
    }
    std::string_view backendId() const noexcept override { return "star-2.7.11b-seed-v1"; }
    const SeedDeliveryContract& outputContract() const noexcept override {
        return plan_->delivery;
    }
    Result<std::unique_ptr<StageWorkspace>> createWorkspace(
        const CompatibilityLimits&, const WorkspaceBudget&) const override {
        return std::unique_ptr<StageWorkspace>(
            std::make_unique<StarSeedWorkspace>(workspace_identity_,
                StarSeedResources{index_, reference_, chromosome_bins_}));
    }
    Status find(const MappingRead& read, StageWorkspace& workspace,
                SeedBatch& out) const override {
        auto* actual = static_cast<StarSeedWorkspace*>(&workspace);
        out.occurrences.reset();
        auto source = actual->sources->acquire(actual->resources, plan_->delivery);
        source->reset(read);
        auto& native = source->native;
        if (auto status = search_.search(read, native); std::holds_alternative<Failure>(status))
            return std::get<Failure>(std::move(status));
        out.seeds.clear();
        out.native_trace.reset();
        out.markers = {};
        out.generation = reference_->view().generation;
        out.splits.assign(native.splits.begin(), native.splits.end());
        out.delivery = plan_->delivery;
        out.markers.legacy_map_marker = native.map_marker;
        out.markers.normal_early_stop = native.map_marker != 0;
        out.occurrences = source;
        out.seeds.reserve(source->pieces().size());
        for (std::size_t i = 0; i < source->pieces().size(); ++i) {
            const auto& piece = source->piece(i);
            out.seeds.push_back(Seed{
                ReadRange{ReadOffset{piece.read_start},
                          ReadOffset{piece.read_start+piece.length}},
                piece.fragment == 0 ? MateId::First : MateId::Second,
                MatchMultiplicity{piece.multiplicity, CountPrecision::Exact},
                OccurrenceHandle{source->identity().source, OccurrenceSetId{i}}});
        }
        return std::monostate{};
    }
private:
    std::shared_ptr<const SuffixArrayIndex> index_;
    std::shared_ptr<const GenomeReferenceData> reference_;
    std::shared_ptr<const SeedSearchPlan> plan_;
    StarSeedSearch search_;
    WorkspaceIdentity workspace_identity_;
    std::shared_ptr<const std::vector<std::uint32_t>> chromosome_bins_;
};

} // namespace

Result<std::shared_ptr<const ISeedFinder>> StarSeedAlgorithmFactory::bind(
    SearchIndexLease index, ReferenceLease reference,
    std::shared_ptr<const SeedSearchPlan> plan) const {
    const auto sa = std::dynamic_pointer_cast<const SuffixArrayIndex>(index);
    const auto genome = std::dynamic_pointer_cast<const GenomeReferenceData>(reference);
    if (!sa || !genome || !plan)
        return seedError(ErrorKind::Config, "STAR seed search needs SA index and reference");
    const auto& delivery = plan->delivery;
    if (delivery.strategy != SeedContract::Star2711bEquivalent ||
        delivery.hit_order != HitOrder::Star2711b ||
        !delivery.exact_multiplicity || !delivery.complete_occurrences ||
        !delivery.known_junction_projection || delivery.hit_schema_version != 1 ||
        delivery.search_domain_id != star_seed_search_domain ||
        !plan->algorithm_options.schema_id.empty() ||
        sa->identity().reference_semantic_id != genome->semanticId())
        return seedError(ErrorKind::Config, "unsupported STAR seed delivery contract");
    try {
        return std::shared_ptr<const ISeedFinder>(
            std::make_shared<StarSeedFinder>(std::move(sa), std::move(genome),
                                             std::move(plan)));
    } catch (const std::exception& error) {
        return seedError(ErrorKind::Config, error.what());
    }
}

} // namespace nustar::star2711b
