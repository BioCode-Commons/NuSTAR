#pragma once

#include "../reference/view.hpp"

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>

namespace nustar::star2711b::detail {

inline std::uint64_t reverseBytes(std::uint64_t word) noexcept {
    word=((word&0x00ff00ff00ff00ffULL)<<8)|((word>>8)&0x00ff00ff00ff00ffULL);
    word=((word&0x0000ffff0000ffffULL)<<16)|((word>>16)&0x0000ffff0000ffffULL);
    return (word<<32)|(word>>32);
}

template<bool Forward>
inline std::uint64_t traversalWord(const BaseCode* bases) noexcept {
    std::uint64_t word;
    std::memcpy(&word,Forward ? bases : bases-7,sizeof(word));
    static_assert(std::endian::native==std::endian::little ||
                  std::endian::native==std::endian::big);
    if constexpr (Forward != (std::endian::native==std::endian::little))
        word=reverseBytes(word);
    return word;
}

// STAR's byte comparisons, grouped into bounded machine-word loads. The first
// differing byte still determines the matched length and ordering, including
// reverse-genome N handling. No read extends beyond the requested query span;
// wide genome loads stay inside the logical span, leaving guards to the tail.
template<bool QueryForward, bool GenomeForward>
inline std::uint64_t compareBases(const BaseCode* query, const BaseCode* genome,
    std::uint64_t length, const BaseCode* genome_begin, const BaseCode* genome_end,
    bool& query_greater) noexcept {
    const auto available=genome>=genome_begin && genome<genome_end ?
        std::uint64_t(GenomeForward ? genome_end-genome : genome-genome_begin+1) : 0;
    const auto wide_end=std::min(length,available)/8*8;
    const auto mismatch=[&](std::uint8_t a,std::uint8_t b) {
        query_greater=GenomeForward ? a>b : !(a>b || b>3);
    };
    std::uint64_t i=0;
    for (;i<wide_end;i+=8) {
        const auto offset=static_cast<std::ptrdiff_t>(i);
        const auto a=traversalWord<QueryForward>(query+(QueryForward ? offset : -offset));
        const auto b=traversalWord<GenomeForward>(genome+(GenomeForward ? offset : -offset));
        if (const auto different=a^b; different) {
            const auto first=std::countr_zero(different)/8;
            mismatch(static_cast<std::uint8_t>(a>>(8*first)),
                     static_cast<std::uint8_t>(b>>(8*first)));
            return i+first;
        }
    }
    for (;i<length;++i) {
        const auto offset=static_cast<std::ptrdiff_t>(i);
        const auto a=static_cast<std::uint8_t>(query[QueryForward ? offset : -offset]);
        const auto b=static_cast<std::uint8_t>(genome[GenomeForward ? offset : -offset]);
        if (a!=b) { mismatch(a,b); return i; }
    }
    return length;
}

} // namespace nustar::star2711b::detail
