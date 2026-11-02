#include "suffix_array_index.hpp"
#include "../reference/genome.hpp"

#include <algorithm>
#include <bit>
#include <limits>
#include <stdexcept>

namespace nustar::star2711b {
namespace {
// Decode the existing physical SA once per distinct prefix, using STAR's
// forward genome and reverse-complement comparison convention.
std::pair<std::uint64_t, unsigned> prefixAt(const PackedSuffixArray& sa,
    std::uint64_t rank, const GenomeReferenceData& reference, unsigned bases) {
    const auto real = reference.storage().forward.size();
    const auto synthetic = reference.synthetic().size();
    const auto genome = reference.starGenome();
    const auto p = sa[rank];
    std::uint64_t start;
    bool reverse;
    if (p < real) { start=p; reverse=false; }
    else if (p >= real+20 && p < 2*real+20) {
        start=2*real+19-p; reverse=true;
    } else if (p >= 2*real+40 && p < 2*real+40+synthetic) {
        start=p-real-40; reverse=false;
    } else if (p >= 2*real+60+synthetic && p < 2*real+60+2*synthetic) {
        start=genome.size()-1-(p-(2*real+60+synthetic)); reverse=true;
    } else throw std::runtime_error("prefix SA points at a separator");
    std::uint64_t code=0;
    unsigned length=0;
    for (; length<bases; ++length) {
        if ((reverse && length>start) || (!reverse && start+length>=genome.size())) break;
        auto base=std::uint8_t(genome[reverse ? start-length : start+length]);
        if (base>3) break;
        if (reverse) base=3-base;
        code=(code<<2)|base;
    }
    return {code,length};
}
}

unsigned StarPrefixIndex::width(std::uint64_t genome_length) {
    // STAR uses at least 32 strand-position bits, followed by rank/N/absent.
    const auto strand=std::max(32U, unsigned(std::bit_width(genome_length)));
    if (strand>61) throw std::length_error("STAR prefix flags exceed 64 bits");
    return strand+3;
}
std::uint64_t StarPrefixIndex::entryCount(std::uint32_t bases) {
    if (!bases || bases>15) throw std::invalid_argument("invalid STAR prefix depth");
    return ((std::uint64_t{1} << (2*(bases+1)))-4)/3;
}
void StarPrefixIndex::layout(unsigned bases, std::uint64_t genome_length) {
    entryCount(bases);
    starts_.assign(bases+1,0);
    for (unsigned i=1;i<=bases;++i) starts_[i]=starts_[i-1]+(std::uint64_t{1}<<(2*i));
    n_mask_=std::uint64_t{1}<<(width(genome_length)-2);
    absent_mask_=n_mask_<<1;
}
StarPrefixIndex::StarPrefixIndex(std::uint32_t bases, std::uint64_t genome_length,
    const PackedSuffixArray& sa, const GenomeReferenceData& reference,
    std::span<const PrefixRange> legacy) {
    layout(bases,genome_length);
    entries_=PackedSuffixArray::zeroed(entryCount(bases),width(genome_length));
    std::vector<std::uint64_t> next(bases,0), genomic_next(bases,0);
    const auto real=reference.storage().forward.size();
    auto genomic=[&](std::uint64_t rank) {
        const auto p=sa[rank];
        return p<real || (p>=real+20 && p<2*real+20);
    };
    auto range=legacy.begin();
    // The legacy table is an acceleration hint for construction: every full
    // ACGT prefix is already a contiguous SA range. Gaps are inspected as in
    // genomeSAindexChunk, so short prefixes and N flags are not discarded.
    for (std::uint64_t rank=0;rank<sa.size();) {
        std::uint64_t code, following;
        unsigned length;
        if (range!=legacy.end() && range->begin==rank) {
            code=range->code; length=bases; following=range->end; ++range;
        } else {
            auto decoded=prefixAt(sa,rank,reference,bases);
            code=decoded.first; length=decoded.second; following=rank+1;
        }
        bool has_genomic=false;
        for (auto i=rank;i<following && !has_genomic;++i) has_genomic=genomic(i);
        for (unsigned level=0;level<length;++level) {
            const auto prefix=code>>(2*(length-1-level));
            if (has_genomic) genomic_next[level]=prefix+1;
            if (prefix>=next[level]) {
                for (auto missing=next[level];missing<prefix;++missing)
                    entries_.set(starts_[level]+missing,rank|absent_mask_);
                entries_.set(starts_[level]+prefix,rank);
                next[level]=prefix+1;
            }
        }
        for (unsigned level=length;level<bases;++level) {
            // STAR marks the preceding valid prefix when suffixes containing
            // N follow it. An empty level has no preceding valid entry.
            // sjdbBuildIndex retains the original genome's N flag even if
            // a newly inserted junction prefix is now closer to that N run.
            if (has_genomic && genomic_next[level]) {
                const auto i=starts_[level]+genomic_next[level]-1;
                entries_.set(i,entries_[i]|n_mask_);
            }
            if (!has_genomic && next[level]) {
                const auto i=starts_[level]+next[level]-1;
                entries_.set(i,entries_[i]|n_mask_);
            }
        }
        rank=following;
    }
    for (unsigned level=0;level<bases;++level)
        for (auto i=starts_[level]+next[level];i<starts_[level+1];++i)
            entries_.set(i,sa.size()|absent_mask_);
}
StarPrefixIndex StarPrefixIndex::fromBytes(std::uint32_t bases,
    std::uint64_t genome_length,
    std::vector<std::uint8_t> bytes) {
    StarPrefixIndex result;
    result.layout(bases,genome_length);
    result.entries_=PackedSuffixArray::fromBytes(width(genome_length),entryCount(bases),std::move(bytes));
    return result;
}
StarPrefixIndex StarPrefixIndex::fromView(std::uint32_t bases,
    std::uint64_t genome_length,
    std::span<const std::uint8_t> bytes, std::shared_ptr<const void> owner) {
    StarPrefixIndex result;
    result.layout(bases,genome_length);
    result.entries_=PackedSuffixArray::fromView(width(genome_length),entryCount(bases),bytes,std::move(owner));
    return result;
}

std::vector<PrefixRange> SuffixArrayIndex::prefixRanges() const {
    // Diagnostics and v2 export only; normal mapping uses the packed SAi view.
    const auto& ref=dynamic_cast<const GenomeReferenceData&>(*reference_);
    const auto begin=prefix_index_.levelStart(prefix_bases_);
    const auto end=prefix_index_.levelEnd(prefix_bases_);
    std::vector<PrefixRange> ranges;
    for (auto i=begin;i<end;++i) {
        const auto entry=prefix_index_[i];
        if (entry&prefix_index_.absentMask()) continue;
        auto first=entry&prefix_index_.rankMask();
        auto last=i+1==end ? suffixes_.size() : prefix_index_[i+1]&prefix_index_.rankMask();
        // Retained pre-insertion N flags need not mark the immediately
        // preceding newly inserted SJ prefix. Check the actual range tail
        // when reconstructing the old diagnostic table; search uses SAi flags.
        const auto tail=last>first ? prefixAt(suffixes_,last-1,ref,prefix_bases_) :
            std::pair<std::uint64_t,unsigned>{0,0};
        if ((entry&prefix_index_.nMask()) || tail.second!=prefix_bases_ || tail.first!=i-begin) {
            // N suffixes sort after the full ACGT prefix; locate its first
            // excluded rank without scanning a large repeated range.
            auto lo=first, hi=last;
            while (lo<hi) {
                const auto mid=lo+(hi-lo)/2;
                const auto p=prefixAt(suffixes_,mid,ref,prefix_bases_);
                if (p.second==prefix_bases_ && p.first==i-begin) lo=mid+1;
                else hi=mid;
            }
            last=lo;
        }
        ranges.push_back({i-begin,first,last});
    }
    return ranges;
}
} // namespace nustar::star2711b
