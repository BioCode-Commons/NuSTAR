#pragma once

#include <cstddef>
#include <cstdint>
#include <istream>
#include <ostream>
#include <string>

namespace nustar::star2711b::binary {

inline void put32(std::ostream& out, std::uint32_t value) {
    for (unsigned i = 0; i < 4; ++i) out.put(char((value >> (8 * i)) & 255));
}
inline void put64(std::ostream& out, std::uint64_t value) {
    for (unsigned i = 0; i < 8; ++i) out.put(char((value >> (8 * i)) & 255));
}
inline bool get32(std::istream& in, std::uint32_t& value) {
    value = 0;
    for (unsigned i = 0; i < 4; ++i) {
        const int byte = in.get();
        if (byte < 0) return false;
        value |= std::uint32_t(byte) << (8 * i);
    }
    return true;
}
inline bool get64(std::istream& in, std::uint64_t& value) {
    value = 0;
    for (unsigned i = 0; i < 8; ++i) {
        const int byte = in.get();
        if (byte < 0) return false;
        value |= std::uint64_t(byte) << (8 * i);
    }
    return true;
}
inline void putString(std::ostream& out, const std::string& value) {
    put32(out, static_cast<std::uint32_t>(value.size()));
    out.write(value.data(), static_cast<std::streamsize>(value.size()));
}
inline bool getString(std::istream& in, std::string& value, std::uint32_t max_size = 1 << 20) {
    std::uint32_t size = 0;
    if (!get32(in, size) || size > max_size) return false;
    value.resize(size);
    return bool(in.read(value.data(), size));
}

} // namespace nustar::star2711b::binary
