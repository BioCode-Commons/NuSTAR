#pragma once

#include "base_compare.hpp"

#if (defined(__x86_64__) || defined(__i386__)) && \
    (defined(__GNUC__) || defined(__clang__))
#include <immintrin.h>
#define NUSTAR_HAS_AVX2_BASE_COMPARE 1
#else
#define NUSTAR_HAS_AVX2_BASE_COMPARE 0
#endif

namespace nustar::star2711b::detail {

// Detect once per process. The caller can select its comparison function when
// binding the seed search, so no CPU feature check is needed in a probe loop.
inline bool avx2Available() noexcept {
#if NUSTAR_HAS_AVX2_BASE_COMPARE
    static const bool available=[] {
        __builtin_cpu_init();
        return bool(__builtin_cpu_supports("avx2"));
    }();
    return available;
#else
    return false;
#endif
}

#if NUSTAR_HAS_AVX2_BASE_COMPARE

[[gnu::target("avx2")]]
inline __m256i reverseBaseVector(__m256i bytes) noexcept {
    // vpshufb reverses each 128-bit lane; swap the lanes for a full reversal.
    const auto order=_mm256_setr_epi8(
        15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,
        15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0);
    const auto reversed_lanes=_mm256_shuffle_epi8(bytes,order);
    return _mm256_permute2x128_si256(reversed_lanes,reversed_lanes,0x01);
}

// Only call on a CPU accepted by avx2Available(). This function
// has its own target attribute: the rest of the binary needs no -mavx2 flag.
// Reads stay within the requested query span AND the logical genomic span.
// The bounded word implementation handles short tails and genome guards.
template<bool QueryForward, bool GenomeForward>
[[gnu::target("avx2")]]
inline std::uint64_t compareBasesAvx2(const BaseCode* query, const BaseCode* genome,
    std::uint64_t length, const BaseCode* genome_begin, const BaseCode* genome_end,
    bool& query_greater) noexcept {
    const auto available=genome>=genome_begin && genome<genome_end ?
        std::uint64_t(GenomeForward ? genome_end-genome : genome-genome_begin+1) : 0;
    const auto wide_end=std::min(length,available)/32*32;
    std::uint64_t i=0;
    for (;i<wide_end;i+=32) {
        const auto offset=static_cast<std::ptrdiff_t>(i);
        const auto* q=query+(QueryForward ? offset : -offset-31);
        const auto* g=genome+(GenomeForward ? offset : -offset-31);
        auto a=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(q));
        auto b=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(g));
        if constexpr (QueryForward != GenomeForward) {
            if constexpr (!QueryForward) a=reverseBaseVector(a);
            else b=reverseBaseVector(b);
        }
        const auto different=~static_cast<std::uint32_t>(
            _mm256_movemask_epi8(_mm256_cmpeq_epi8(a,b)));
        if (different) {
            // Two backwards operands already share byte order; read their mask
            // from its high bit instead of shuffling both vectors.
            const auto first=(!QueryForward && !GenomeForward) ?
                std::countl_zero(different) : std::countr_zero(different);
            const auto position=static_cast<std::ptrdiff_t>(i+first);
            const auto x=std::uint8_t(query[QueryForward ? position : -position]);
            const auto y=std::uint8_t(genome[GenomeForward ? position : -position]);
            query_greater=GenomeForward ? x>y : !(x>y || y>3);
            return i+first;
        }
    }
    if (i==length) return length;
    const auto offset=static_cast<std::ptrdiff_t>(i);
    return i+compareBases<QueryForward,GenomeForward>(
        query+(QueryForward ? offset : -offset),
        genome+(GenomeForward ? offset : -offset),length-i,
        genome_begin,genome_end,query_greater);
}

#else

// Keep the interface available for portable builds; capability detection is
// false here, and the implementation uses the existing bounded word kernel.
template<bool QueryForward, bool GenomeForward>
inline std::uint64_t compareBasesAvx2(const BaseCode* query, const BaseCode* genome,
    std::uint64_t length, const BaseCode* genome_begin, const BaseCode* genome_end,
    bool& query_greater) noexcept {
    return compareBases<QueryForward,GenomeForward>(query,genome,length,
                                                    genome_begin,genome_end,query_greater);
}

#endif

} // namespace nustar::star2711b::detail
