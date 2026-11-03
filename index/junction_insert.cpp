#include "junction_insert.hpp"
#include <algorithm>
#include <array>
#include <limits>
#include <tuple>

namespace nustar::star2711b {
namespace {
using U = std::uint64_t;
struct Insertion { U rank, position; };
struct Location { U position; bool reverse; };
Location locate(U p, U real, U sj) {
    if (p < real) return {p, false};
    if (p < 2*real+20) return {2*real+19-p, true};
    if (p < 2*real+40+sj) return {p-real-40, false};
    return {3*real+2*sj+59-p, true};
}
// compareSeqToGenome1 / suffixArraySearch1, with gInsert=UINT64_MAX.
// At equal padding, forward old suffixes precede the new suffix, reverse
// old suffixes follow it. An ordinary lexical upper_bound is not equivalent.
U insertionRank(const std::vector<std::uint8_t>& sj, U offset,
    const PackedSuffixArray& sa, const GenomeReferenceData& reference) {
    const U real=reference.storage().forward.size(), syn=reference.synthetic().size();
    const auto genome=reference.starGenome();
    auto compare=[&](U rank, U matched) {
        const auto location=locate(sa[rank],real,syn);
        for (;matched<10000;++matched) {
            const auto q=sj[offset+matched];
            auto g=static_cast<std::uint8_t>(genome[location.reverse ?
                location.position-matched : location.position+matched]);
            if (location.reverse && g<4) g=3-g;
            if (q!=g) return std::pair{matched,q>g ? 1 : -1};
            if (q==5) return std::pair{matched,location.reverse ? -1 : 1};
        }
        return std::pair{matched,0};
    };
    if (sa.empty()) return 0;
    U left=0,right=sa.size()-1;
    auto [ll,lc]=compare(left,0);
    if (lc<0) return 0;
    auto [rl,rc]=compare(right,0);
    if (rc>0) return U(-2);
    while (left+1<right) {
        const U middle=(left+right)/2;
        auto [length,order]=compare(middle,std::min(ll,rl));
        if (length==10000) return middle;
        if (order>0) { left=middle; ll=length; }
        else { right=middle; rl=length; }
    }
    return right;
}
}
Result<SearchIndexLease> insertJunctions(const SuffixArrayIndex& old_index,
    const GenomeReferenceData& old_reference,
    std::shared_ptr<const GenomeReferenceData> reference, unsigned threads) {
    try {
        const U real=reference->storage().forward.size();
        const U stride=2*reference->sjdbOverhang()+1, size=reference->synthetic().size();
        const auto entries=reference->junctions(), previous=old_reference.junctions();
        const auto& old_sa=old_index.suffixArray();
        std::vector<std::uint8_t> sequence(2*size+1,5);
        for(U i=0;i<size;++i) {
            auto b=static_cast<std::uint8_t>(reference->synthetic()[i]);
            sequence[i]=b; sequence[2*size-1-i]=b<4 ? 3-b : b;
        }
        std::vector<U> old_to_new(previous.size());
        std::vector<bool> existing(entries.size(),false);
        for(U i=0;i<entries.size();++i) {
            const auto& j=entries[i];
            auto it=std::lower_bound(previous.begin(),previous.end(),j,[](auto& a,auto& b) {
                return std::tie(a.lookup_intron.begin.value,a.lookup_intron.end.value)<
                       std::tie(b.lookup_intron.begin.value,b.lookup_intron.end.value);
            });
            if(it!=previous.end() && it->lookup_intron.begin.value==j.lookup_intron.begin.value &&
               it->lookup_intron.end.value==j.lookup_intron.end.value) {
                old_to_new[it-previous.begin()]=i; existing[i]=true;
            }
        }
        // Only novel SJ suffixes need insertion records. Existing junctions
        // still use old_to_new, but do not allocate a placeholder per base.
        std::vector<U> tasks;
        for(U block=0;block<2*entries.size();++block) {
            const U logical=block<entries.size() ? block : 2*entries.size()-1-block;
            if(!existing[logical]) tasks.push_back(block);
        }
        std::vector<Insertion> insertions(tasks.size()*stride,{U(-1),0});
        #pragma omp parallel for num_threads(threads) schedule(dynamic,1)
        for(U task=0;task<tasks.size();++task) {
            const U block=tasks[task];
            for(U i=0;i<stride;++i) {
                const U p=block*stride+i;
                if(sequence[p]<4) insertions[task*stride+i]={insertionRank(sequence,p,old_sa,old_reference),p};
            }
        }
        insertions.erase(std::remove_if(insertions.begin(),insertions.end(),
            [](auto i){return i.rank==U(-1);}),insertions.end());
        std::sort(insertions.begin(),insertions.end(),[&](auto a,auto b){
            if(a.rank!=b.rank) return a.rank<b.rank;
            U k=0;
            while(sequence[a.position+k]==sequence[b.position+k]) {
                if(sequence[a.position+k]==5) return a.position<b.position;
                ++k;
            }
            return sequence[a.position+k]<sequence[b.position+k];
        });
        const U physical=2*(real+size)+80;
        auto sa=PackedSuffixArray::lazyZeroed(old_sa.size()+insertions.size(),PackedSuffixArray::widthFor(physical));
        auto newPosition=[&](U p) {return p<size ? 2*real+40+p : 2*real+60+p;};
        const auto old_size=old_reference.synthetic().size();
        auto oldPosition=[&](U p) {
            if(p>=2*real+40 && p<2*real+40+old_size) {
                const U relative=p-(2*real+40);
                p=2*real+40+old_to_new[relative/stride]*stride+relative%stride;
            } else if(p>=2*real+60+old_size) {
                const U relative=p-(2*real+60+old_size);
                const U old_id=previous.size()-1-relative/stride;
                p=2*real+60+size+(entries.size()-1-old_to_new[old_id])*stride+relative%stride;
            }
            return p;
        };
        // Old suffixes keep their order and each insertion goes before the old
        // suffix of its rank (after all of them for U(-2)), so insertion k lands
        // at output rank min(rank,old_count)+k. That rank grows with k, which
        // lets each chunk find its first insertion by binary search and merge
        // its output ranks without the chunks before it.
        const U old_count=old_sa.size(), total=old_count+insertions.size();
        auto insertionsBefore=[&](U output_rank) {
            U low=0,high=insertions.size();
            while(low<high) {
                const U k=low+(high-low)/2;
                if(std::min(insertions[k].rank,old_count)+k<output_rank) low=k+1; else high=k;
            }
            return low;
        };
        // A multiple of 64 ranks, so the chunks' RangeWriters share no byte.
        constexpr U chunk=U(1)<<22;
        const U chunks=(total+chunk-1)/chunk;
        #pragma omp parallel for num_threads(threads) schedule(dynamic,1)
        for(U c=0;c<chunks;++c) {
            const U begin=c*chunk,end=std::min(total,begin+chunk);
            U next=insertionsBefore(begin),rank=begin-next;
            const U first_old=rank;
            auto writer=sa.rangeWriter(begin);
            for(U out=begin;out<end;++out) {
                if(next<insertions.size() && std::min(insertions[next].rank,old_count)==rank)
                    writer.put(newPosition(insertions[next++].position));
                else
                    writer.put(oldPosition(old_sa[rank++]));
            }
            writer.finish();
            old_sa.discardReadPages(first_old,rank);
        }
        // Chunk discards keep pages shared by neighbouring chunks; drop them now.
        old_sa.discardReadPages(0,old_count);

        // Port of the SAi update in sjdbBuildIndex; preserves old N flags.
        const auto& old_prefix=old_index.prefixIndex();
        const auto depth=old_index.prefixBases();
        const auto width=StarPrefixIndex::width(real+size);
        if(width!=StarPrefixIndex::width(real+old_reference.synthetic().size()))
            throw std::runtime_error("cannot insert junctions across STAR strand-bit boundary");
        // Reserve the guards up front: growing a human SAi vector by seven
        // bytes would otherwise double its capacity and briefly keep two copies.
        std::vector<std::uint8_t> prefix_bytes(old_prefix.bytes().size()+7);
        std::copy(old_prefix.bytes().begin(),old_prefix.bytes().end(),prefix_bytes.begin());
        prefix_bytes.resize(old_prefix.bytes().size());
        auto prefix=PackedSuffixArray::fromBytes(width,StarPrefixIndex::entryCount(depth),std::move(prefix_bytes));
        const U absent=old_prefix.absentMask(), nflag=old_prefix.nMask();
        auto codeAt=[&](U pos,unsigned level) -> std::int64_t {
            U code=0;
            for(unsigned k=0;k<=level;++k) {
                if(sequence[pos+k]>3) return -static_cast<std::int64_t>(code);
                code=(code<<2)|sequence[pos+k];
            }
            return code;
        };
        const U count=insertions.size();
        for(unsigned level=0;level<depth;++level) {
            U sj=0, present=old_prefix.levelStart(level+1)-1;
            const U begin=old_prefix.levelStart(level+1),end=old_prefix.levelEnd(level+1);
            for(U i=begin;i<end;++i) {
                const U v=prefix[i],rank=v&old_prefix.rankMask();
                if(sj<count && (v&absent)) {
                    const U saved=sj;
                    auto code=codeAt(insertions[sj].position,level);
                    while(sj<count && code<static_cast<std::int64_t>(i-begin) && insertions[sj].rank-1<rank) {
                        ++sj;
                        if(sj<count) code=codeAt(insertions[sj].position,level);
                    }
                    if(sj<count && code==static_cast<std::int64_t>(i-begin)) {
                        const U value=insertions[sj].rank+sj;
                        prefix.set(i,value);
                        for(U k=present+1;k<i;++k) prefix.set(k,value|absent);
                        ++sj; present=i;
                    } else sj=saved;
                } else {
                    while(sj<count && insertions[sj].rank<rank) ++sj;
                    while(sj<count && insertions[sj].rank==rank) {
                        if(codeAt(insertions[sj].position,level)>=static_cast<std::int64_t>(i-begin)) break;
                        ++sj;
                    }
                    prefix.set(i,v+sj);
                    for(U k=present+1;k<i;++k) prefix.set(k,(rank+sj)|absent);
                    present=i;
                }
            }
        }
        for(const auto insertion:insertions) {
            U code=0;
            for(unsigned level=0;level<depth;++level) {
                const auto base=sequence[insertion.position+level];
                code<<=2;
                if(base>3) {
                    for(unsigned l=level;l<depth;++l) {
                        code+=3;
                        auto i=static_cast<std::int64_t>(old_prefix.levelStart(l+1)+code);
                        while(i>=0 && (prefix[i]&absent)) --i;
                        if(i>=0) prefix.set(i,prefix[i]|nflag);
                        code<<=2;
                    }
                    break;
                }
                code+=base;
            }
        }
        auto owner=std::make_shared<PackedSuffixArray>(std::move(prefix));
        auto native=StarPrefixIndex::fromView(depth,real+size,std::span<const std::uint8_t>(owner->packedBytes().data(),owner->storageBytes()),owner);
        auto identity=old_index.identity();
        identity.reference_semantic_id=reference->view().semantic_id;
        identity.artifact_sha256.clear();
        ++identity.generation.value;
        return SearchIndexLease{std::make_shared<SuffixArrayIndex>(std::move(identity),reference,
            std::filesystem::path{},std::string(old_index.algorithmId()),depth,old_index.suffixLengthMax(),
            physical,std::move(sa),std::vector<PrefixRange>{},std::move(native))};
    } catch(const std::exception& e) {return Failure{ErrorKind::Index,Stage::Sjdb,{},1,e.what()};}
}
}
