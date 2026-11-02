#include "sha256.hpp"

#include <algorithm>
#include <iomanip>
#include <sstream>

#if (defined(__x86_64__) || defined(__i386__)) && \
    (defined(__GNUC__) || defined(__clang__)) && !defined(NUSTAR_SHA256_PORTABLE)
#define NUSTAR_SHA256_X86 1
#include <immintrin.h>
#endif

namespace nustar::star2711b {
namespace {

constexpr std::array<std::uint32_t, 64> k{
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
};

std::uint32_t rotate(std::uint32_t value, unsigned count) {
    return (value >> count) | (value << (32 - count));
}

#ifdef NUSTAR_SHA256_X86
// SHA-NI implements the same SHA-256 rounds and schedule as block(). The
// target is confined to this function; the binary still runs on older CPUs.
// Instruction semantics: Intel SHA Extensions programming reference.
__attribute__((target("sha,ssse3,sse4.1")))
void hardwareBlocks(std::array<std::uint32_t, 8>& state,
                    const std::byte* bytes, std::size_t count) {
    auto abef = _mm_set_epi32(state[0], state[1], state[4], state[5]);
    auto cdgh = _mm_set_epi32(state[2], state[3], state[6], state[7]);
    const auto endian = _mm_set_epi8(12,13,14,15,8,9,10,11,4,5,6,7,0,1,2,3);
    for (std::size_t block = 0; block < count; ++block, bytes += 64) {
        const auto saved_abef = abef, saved_cdgh = cdgh;
        __m128i words[4];
        for (unsigned group = 0; group < 16; ++group) {
            const auto slot = group%4;
            if (group < 4) {
                words[slot] = _mm_shuffle_epi8(
                    _mm_loadu_si128(reinterpret_cast<const __m128i*>(bytes+16*group)), endian);
            } else {
                auto next = _mm_sha256msg1_epu32(words[slot], words[(slot+1)%4]);
                next = _mm_add_epi32(next,
                    _mm_alignr_epi8(words[(slot+3)%4], words[(slot+2)%4], 4));
                words[slot] = _mm_sha256msg2_epu32(next, words[(slot+3)%4]);
            }
            auto message = _mm_add_epi32(words[slot],
                _mm_loadu_si128(reinterpret_cast<const __m128i*>(k.data()+group*4)));
            cdgh = _mm_sha256rnds2_epu32(cdgh, abef, message);
            message = _mm_shuffle_epi32(message, 0x0e);
            abef = _mm_sha256rnds2_epu32(abef, cdgh, message);
        }
        abef = _mm_add_epi32(abef, saved_abef);
        cdgh = _mm_add_epi32(cdgh, saved_cdgh);
    }
    std::uint32_t a[4], b[4];
    _mm_storeu_si128(reinterpret_cast<__m128i*>(a), abef);
    _mm_storeu_si128(reinterpret_cast<__m128i*>(b), cdgh);
    state = {a[3], a[2], b[3], b[2], a[1], a[0], b[1], b[0]};
}
#endif

} // namespace

Sha256::Sha256()
    : state_{0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
             0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19} {}

void Sha256::block(const std::byte* bytes) {
    std::array<std::uint32_t, 64> words{};
    for (std::size_t i = 0; i < 16; ++i) {
        for (std::size_t j = 0; j < 4; ++j) {
            words[i] = (words[i] << 8) | std::to_integer<std::uint8_t>(bytes[4 * i + j]);
        }
    }
    for (std::size_t i = 16; i < words.size(); ++i) {
        const auto a = rotate(words[i - 15], 7) ^ rotate(words[i - 15], 18) ^ (words[i - 15] >> 3);
        const auto b = rotate(words[i - 2], 17) ^ rotate(words[i - 2], 19) ^ (words[i - 2] >> 10);
        words[i] = words[i - 16] + a + words[i - 7] + b;
    }
    auto [a,b,c,d,e,f,g,h] = state_;
    for (std::size_t i = 0; i < words.size(); ++i) {
        const auto s1 = rotate(e,6) ^ rotate(e,11) ^ rotate(e,25);
        const auto choice = (e & f) ^ (~e & g);
        const auto t1 = h + s1 + choice + k[i] + words[i];
        const auto s0 = rotate(a,2) ^ rotate(a,13) ^ rotate(a,22);
        const auto majority = (a & b) ^ (a & c) ^ (b & c);
        const auto t2 = s0 + majority;
        h=g; g=f; f=e; e=d+t1; d=c; c=b; b=a; a=t1+t2;
    }
    state_[0]+=a; state_[1]+=b; state_[2]+=c; state_[3]+=d;
    state_[4]+=e; state_[5]+=f; state_[6]+=g; state_[7]+=h;
}

void Sha256::update(std::span<const std::byte> bytes) {
    bytes_ += bytes.size();
    while (!bytes.empty()) {
        if (buffered_ == 0 && bytes.size() >= buffer_.size()) {
            const auto count = bytes.size()/buffer_.size();
            blocks(bytes.data(), count);
            bytes = bytes.subspan(count*buffer_.size());
            continue;
        }
        const auto count = std::min(bytes.size(), buffer_.size() - buffered_);
        std::copy_n(bytes.begin(), count, buffer_.begin() + buffered_);
        bytes = bytes.subspan(count);
        buffered_ += count;
        if (buffered_ == buffer_.size()) {
            blocks(buffer_.data(), 1);
            buffered_ = 0;
        }
    }
}

void Sha256::blocks(const std::byte* bytes, std::size_t count) {
#ifdef NUSTAR_SHA256_X86
    static const bool accelerated = __builtin_cpu_supports("sha") &&
        __builtin_cpu_supports("ssse3") && __builtin_cpu_supports("sse4.1");
    if (accelerated) {
        hardwareBlocks(state_, bytes, count);
        return;
    }
#endif
    for (std::size_t i = 0; i < count; ++i) block(bytes+64*i);
}

void Sha256::update(std::string_view text) {
    update(std::as_bytes(std::span(text.data(), text.size())));
}

std::string Sha256::finish() {
    const auto bit_count = bytes_ * 8;
    const std::byte one{0x80};
    update(std::span(&one, 1));
    const std::byte zero{0};
    while (buffered_ != 56) update(std::span(&zero, 1));
    std::array<std::byte, 8> length{};
    for (std::size_t i = 0; i < length.size(); ++i) {
        length[7 - i] = std::byte((bit_count >> (8 * i)) & 0xff);
    }
    update(length);
    std::ostringstream out;
    out << std::hex << std::setfill('0');
    for (const auto word : state_) out << std::setw(8) << word;
    return out.str();
}

std::string sha256(std::span<const std::byte> bytes) {
    Sha256 hash;
    hash.update(bytes);
    return hash.finish();
}

} // namespace nustar::star2711b
