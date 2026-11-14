#pragma once

#include "base_compare.hpp"

#if defined(__aarch64__) && defined(__ARM_NEON)
#include <arm_neon.h>
#endif

namespace nustar::star2711b::detail {

#if defined(__aarch64__) && defined(__ARM_NEON)
inline constexpr bool neonBaseComparisonAvailable = true;

inline uint8x16_t reverseNeon(uint8x16_t value) noexcept {
    // Reverse bytes within each half, then swap the two halves.
    value = vrev64q_u8(value);
    return vextq_u8(value, value, 8);
}
#else
inline constexpr bool neonBaseComparisonAvailable = false;
#endif

// AArch64 NEON counterpart of compareBases. Full 16-byte loads are restricted
// to both the requested query length and the logical genome span; the existing
// bounded word/scalar implementation handles tails and genome guard bytes.
// No complement is applied here: the caller supplies the appropriate encoded
// query, just as it does for compareBases.
template<bool QueryForward, bool GenomeForward>
inline std::uint64_t compareBasesNeon(const BaseCode* query, const BaseCode* genome,
    std::uint64_t length, const BaseCode* genome_begin, const BaseCode* genome_end,
    bool& query_greater) noexcept {
#if defined(__aarch64__) && defined(__ARM_NEON)
    const auto available = genome >= genome_begin && genome < genome_end ?
        std::uint64_t(GenomeForward ? genome_end-genome : genome-genome_begin+1) : 0;
    const auto wide_end = std::min(length, available)/16*16;
    std::uint64_t i = 0;
    for (; i < wide_end; i += 16) {
        const auto offset = static_cast<std::ptrdiff_t>(i);
        const auto* q = query+(QueryForward ? offset : -offset-15);
        const auto* g = genome+(GenomeForward ? offset : -offset-15);
        auto a = vld1q_u8(reinterpret_cast<const std::uint8_t*>(q));
        auto b = vld1q_u8(reinterpret_cast<const std::uint8_t*>(g));
        if constexpr (QueryForward != GenomeForward) {
            if constexpr (!QueryForward) a = reverseNeon(a);
            else b = reverseNeon(b);
        }
        const auto different = veorq_u8(a, b);
        if (vmaxvq_u8(different)) {
            auto low = vgetq_lane_u64(vreinterpretq_u64_u8(different), 0);
            auto high = vgetq_lane_u64(vreinterpretq_u64_u8(different), 1);
            if constexpr (std::endian::native == std::endian::big) {
                low = reverseBytes(low);
                high = reverseBytes(high);
            }
            // Two backwards operands already share byte order. Start from
            // the high byte instead of reversing both vectors on every load.
            unsigned first;
            if constexpr (!QueryForward && !GenomeForward)
                first = high ? std::countl_zero(high)/8 : 8+std::countl_zero(low)/8;
            else
                first = low ? std::countr_zero(low)/8 : 8+std::countr_zero(high)/8;
            const auto mismatch_offset = static_cast<std::ptrdiff_t>(i+first);
            const auto q = static_cast<std::uint8_t>(
                query[QueryForward ? mismatch_offset : -mismatch_offset]);
            const auto g = static_cast<std::uint8_t>(
                genome[GenomeForward ? mismatch_offset : -mismatch_offset]);
            query_greater = GenomeForward ? q > g : !(q > g || g > 3);
            return i+first;
        }
    }
    if (i == length) return length; // Preserve query_greater on exact equality.
    const auto offset = static_cast<std::ptrdiff_t>(i);
    return i+compareBases<QueryForward, GenomeForward>(
        query+(QueryForward ? offset : -offset),
        genome+(GenomeForward ? offset : -offset), length-i,
        genome_begin, genome_end, query_greater);
#else
    return compareBases<QueryForward, GenomeForward>(query, genome, length,
        genome_begin, genome_end, query_greater);
#endif
}

} // namespace nustar::star2711b::detail
