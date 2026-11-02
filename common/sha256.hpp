#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <string_view>

namespace nustar::star2711b {

class Sha256 {
public:
    Sha256();
    void update(std::span<const std::byte> bytes);
    void update(std::string_view text);
    std::string finish();
private:
    void block(const std::byte* bytes);
    void blocks(const std::byte* bytes, std::size_t count);
    std::array<std::uint32_t, 8> state_;
    std::array<std::byte, 64> buffer_{};
    std::uint64_t bytes_ = 0;
    std::size_t buffered_ = 0;
};

std::string sha256(std::span<const std::byte> bytes);

} // namespace nustar::star2711b
