#pragma once

#include "base_compare_avx2.hpp"
#include "base_compare_neon.hpp"

namespace nustar::star2711b::detail {

// Most interval probes fail immediately. Check the first two bases before
// paying for vector setup; matching probes use bounded SIMD/word loads.
template<bool QueryForward, bool GenomeForward>
inline std::uint64_t compareSeedBases(const BaseCode* query, const BaseCode* genome,
    std::uint64_t length, const BaseCode* genome_begin, const BaseCode* genome_end,
    bool& query_greater, bool use_avx2) noexcept {
    const auto head=std::min<std::uint64_t>(2,length);
    for (std::uint64_t i=0;i<head;++i) {
        const auto offset=static_cast<std::ptrdiff_t>(i);
        const auto a=std::uint8_t(query[QueryForward ? offset : -offset]);
        const auto b=std::uint8_t(genome[GenomeForward ? offset : -offset]);
        if (a!=b) {
            query_greater=GenomeForward ? a>b : !(a>b || b>3);
            return i;
        }
    }
    if (head==length) return length;
    const auto offset=static_cast<std::ptrdiff_t>(head);
    query+=QueryForward ? offset : -offset;
    genome+=GenomeForward ? offset : -offset;
    if constexpr (neonBaseComparisonAvailable)
        return head+compareBasesNeon<QueryForward,GenomeForward>(
            query,genome,length-head,genome_begin,genome_end,query_greater);
    if (use_avx2 && length-head>=32)
        return head+compareBasesAvx2<QueryForward,GenomeForward>(
            query,genome,length-head,genome_begin,genome_end,query_greater);
    return head+compareBases<QueryForward,GenomeForward>(
        query,genome,length-head,genome_begin,genome_end,query_greater);
}

} // namespace nustar::star2711b::detail
