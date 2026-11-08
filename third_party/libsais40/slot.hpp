#pragma once
#include <bit>
#include <cstdint>
#include <cstring>
#include <type_traits>

namespace nustar::star2711b::detail::libsais40_impl {
constexpr std::uint64_t slot_mask = (std::uint64_t{1} << 40) - 1;
constexpr std::int64_t signed40(std::uint64_t value) noexcept {
    value &= slot_mask;
    return std::bit_cast<std::int64_t>(value | ((value & (std::uint64_t{1} << 39)) ? ~slot_mask : 0));
}
constexpr std::uint64_t unsigned40(std::uint64_t value) noexcept { return value & slot_mask; }
// Only array entries are packed. Arithmetic, lengths and bounded caches stay native.
struct Slot {
    std::uint8_t bytes[5];
    operator std::int64_t() const noexcept {
        std::uint32_t low;
        std::memcpy(&low, bytes, 4);
        if constexpr (std::endian::native == std::endian::big) low = __builtin_bswap32(low);
        return signed40(low | (std::uint64_t(bytes[4]) << 32));
    }
    Slot& operator=(std::int64_t value) noexcept {
        auto low = static_cast<std::uint32_t>(value);
        if constexpr (std::endian::native == std::endian::big) low = __builtin_bswap32(low);
        std::memcpy(bytes, &low, 4);
        bytes[4] = static_cast<std::uint64_t>(value) >> 32;
        return *this;
    }
    std::int64_t operator++(int) noexcept { auto v = std::int64_t(*this); *this = v + 1; return v; }
    std::int64_t operator--(int) noexcept { auto v = std::int64_t(*this); *this = v - 1; return v; }
    Slot& operator++() noexcept { return *this = std::int64_t(*this) + 1; }
    Slot& operator--() noexcept { return *this = std::int64_t(*this) - 1; }
    Slot& operator+=(std::int64_t v) noexcept { return *this = std::int64_t(*this) + v; }
    Slot& operator-=(std::int64_t v) noexcept { return *this = std::int64_t(*this) - v; }
    Slot& operator|=(std::int64_t v) noexcept { return *this = std::int64_t(*this) | v; }
    Slot& operator&=(std::int64_t v) noexcept { return *this = std::int64_t(*this) & v; }
    Slot& operator^=(std::int64_t v) noexcept { return *this = std::int64_t(*this) ^ v; }
};
static_assert(sizeof(Slot) == 5 && alignof(Slot) == 1);
static_assert(std::is_trivially_copyable_v<Slot> && std::is_trivially_default_constructible_v<Slot>);

// Three-slot scratch records remain Slot objects throughout. A view avoids
// aliasing the SA as a differently aligned native C struct.
struct Record { std::int64_t name,cursor,index; };
struct RecordRef {
    Slot& name;
    Slot& cursor;
    Slot& index;
    operator Record() const noexcept { return {name,cursor,index}; }
    RecordRef& operator=(Record value) noexcept {
        name=value.name; cursor=value.cursor; index=value.index; return *this;
    }
    RecordRef& operator=(const RecordRef& value) noexcept {
        return *this=static_cast<Record>(value);
    }
};
struct RecordView {
    Slot* slots;
    RecordRef operator[](std::ptrdiff_t i) const noexcept {
        return {slots[3*i],slots[3*i+1],slots[3*i+2]};
    }
    RecordView operator+(std::ptrdiff_t i) const noexcept { return {slots+3*i}; }
};
}
