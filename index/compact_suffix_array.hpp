#pragma once
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <span>
#include <memory>
#include <utility>

namespace nustar::star2711b::detail {
// Byte-addressable slots: adjacent writers never modify overlapping bytes.
struct Position40 {
    std::uint8_t bytes[5];
    static constexpr std::uint64_t flag = std::uint64_t{1} << 39;
    static constexpr std::uint64_t mask = 2*flag-1;
    operator std::uint64_t() const noexcept {
        std::uint32_t low; std::memcpy(&low,bytes,4);
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
        low=__builtin_bswap32(low);
#endif
        return low | (std::uint64_t(bytes[4])<<32);
    }
    Position40& operator=(std::uint64_t value) noexcept {
        assert(value<=mask);
        auto low=static_cast<std::uint32_t>(value);
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
        low=__builtin_bswap32(low);
#endif
        std::memcpy(bytes,&low,4); bytes[4]=value>>32; return *this;
    }
};
static_assert(sizeof(Position40)==5 && alignof(Position40)==1);

class CompactSuffixArray {
public:
    using value_type=Position40;
    explicit CompactSuffixArray(std::size_t count);
    ~CompactSuffixArray();
    CompactSuffixArray(CompactSuffixArray&&) noexcept;
    CompactSuffixArray& operator=(CompactSuffixArray&&) noexcept;
    CompactSuffixArray(const CompactSuffixArray&)=delete;
    CompactSuffixArray& operator=(const CompactSuffixArray&)=delete;
    std::span<Position40> slots() noexcept { return {data_,count_}; }
    std::span<const Position40> slots() const noexcept { return {data_,count_}; }
    std::size_t size() const noexcept { return count_; }
    std::size_t allocatedBytes() const noexcept { return bytes_; }
    Position40* begin() noexcept { return data_; }
    Position40* end() noexcept { return data_+count_; }
    Position40& operator[](std::size_t i) noexcept { return data_[i]; }
    std::uint64_t operator[](std::size_t i) const noexcept { return data_[i]; }
    void resize(std::size_t n) noexcept { assert(n<=count_); count_=n; }
    // Only fully consumed pages are discarded; storage must not be read again.
    void discardBefore(std::size_t count) noexcept;
    // Transfer the allocation after in-place packing; unmap unused tail pages.
    // retained_bytes includes the packed reader's seven guard bytes.
    std::shared_ptr<const void> takePackedStorage(std::size_t retained_bytes);
private:
    void release() noexcept;
    Position40* data_=nullptr;
    std::size_t count_=0,bytes_=0,released_=0;
};
CompactSuffixArray buildSuffixArrayLibsais40(std::span<const std::uint8_t>, int runThreadN=0);
void validateLibsais40Length(std::size_t length);
void normalizeStarSuffixRanks40(std::span<Position40>, std::span<const std::uint8_t>,
    std::uint64_t real, std::uint64_t synthetic, std::uint32_t limit,int runThreadN=0);
} // namespace
