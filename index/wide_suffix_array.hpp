#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace nustar::star2711b::detail {

// Internal construction helpers. Search/index interfaces and persisted v3
// packing are independent of this temporary, full-width sorting workspace.
int suffixConstructionThreads(std::size_t length, int runThreadN = 0);
std::vector<std::uint64_t> buildSuffixArraySais64(
    std::span<const std::uint8_t> symbols, int runThreadN = 0);
void normalizeStarSuffixRanks64(std::span<std::uint64_t> suffixes,
    std::span<const std::uint8_t> sequence, std::uint64_t real_length,
    std::uint64_t synthetic_length, std::uint32_t suffix_length_max, int runThreadN = 0);

} // namespace nustar::star2711b::detail
