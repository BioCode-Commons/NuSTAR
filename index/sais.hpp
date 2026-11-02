#pragma once

#include "suffix_array.hpp"

namespace nustar::star2711b {

// Symbols must be positive. A unique zero sentinel is appended internally.
std::vector<std::uint32_t> buildSuffixArraySais(
    std::span<const std::uint8_t> symbols);

class SaisSuffixArrayBuilder final : public ISuffixArrayBuilder {
public:
    std::string_view algorithmId() const noexcept override { return "sais-v1"; }
    std::vector<std::uint32_t> build(
        std::span<const std::uint8_t> symbols) const override {
        return buildSuffixArraySais(symbols);
    }
};

} // namespace nustar::star2711b
