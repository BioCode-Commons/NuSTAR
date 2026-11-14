#pragma once

#include "../reference/view.hpp"

#include <bit>
#include <cstdint>
#include <cstring>

namespace nustar::star2711b::detail {

inline std::uint64_t scalarPrefixCode(const BaseCode* query,unsigned length) noexcept {
    std::uint64_t code=0;
    for (unsigned i=0;i<length;++i) code=(code<<2)|std::uint8_t(query[i]);
    return code;
}

inline std::uint32_t prefixLittleWord(const BaseCode* query) noexcept {
    std::uint32_t word;
    std::memcpy(&word,query,sizeof(word));
    static_assert(std::endian::native==std::endian::little ||
                  std::endian::native==std::endian::big);
    if constexpr (std::endian::native==std::endian::big)
        word=((word&0xffU)<<24)|((word&0xff00U)<<8)|
             ((word>>8)&0xff00U)|(word>>24);
    return word;
}

// Four canonical base bytes b0..b3 become b0<<6 | b1<<4 | b2<<2 | b3.
// In the product, the desired terms occupy disjoint bits [32,39]. All other
// terms start at bit >=42 or have a combined maximum below 2^30, so neither
// overlap nor carry can change those eight bits. The caller checks each byte
// is in [0,3].
inline std::uint64_t packFourPrefixBases(std::uint32_t word) noexcept {
    return (std::uint64_t(word)*0x4010040100ULL>>32)&0xff;
}

// A bounded, allocation-free encoder for SAi prefixes (normally <=15 bases).
// Each load covers exactly four available bytes. Noncanonical base codes take
// the original shift/OR path, preserving its N/spacer semantics exactly.
inline std::uint64_t prefixCode(const BaseCode* query,unsigned length) noexcept {
    if (length<4 || length>15) return scalarPrefixCode(query,length);
    const auto a=prefixLittleWord(query);
    std::uint32_t b=0,c=0;
    if (length>=8) b=prefixLittleWord(query+4);
    if (length>=12) c=prefixLittleWord(query+8);
    if ((a|b|c)&0xfcfcfcfcU) return scalarPrefixCode(query,length);
    auto code=packFourPrefixBases(a);
    unsigned consumed=4;
    if (length>=8) { code=(code<<8)|packFourPrefixBases(b);consumed=8; }
    if (length>=12) { code=(code<<8)|packFourPrefixBases(c);consumed=12; }
    for (;consumed<length;++consumed)
        code=(code<<2)|std::uint8_t(query[consumed]);
    return code;
}

} // namespace nustar::star2711b::detail
