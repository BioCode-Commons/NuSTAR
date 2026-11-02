#pragma once

#include <cstdint>
#include <memory>
#include <span>
#include <string_view>
#include <vector>

namespace nustar::star2711b {

// Builds a lexicographically sorted full suffix array, including the terminal
// suffix at symbols.size(). algorithmId() must remain stable for persisted files.
// The index backend owns layout, STAR rank normalization, prefix ranges and files.
class ISuffixArrayBuilder {
public:
    virtual ~ISuffixArrayBuilder() = default;
    virtual std::string_view algorithmId() const noexcept = 0;
    virtual std::vector<std::uint32_t> build(
        std::span<const std::uint8_t> symbols) const = 0;
};
using SuffixArrayBuilderLease = std::shared_ptr<const ISuffixArrayBuilder>;

} // namespace nustar::star2711b
