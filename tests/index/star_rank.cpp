#include "../../index/wide_suffix_array.hpp"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

int main() {
    // Long exact repeats and long common prefixes ending in a mismatch. The
    // oracle sorts explicit terminated strings, independently of SA/LCP logic.
    constexpr std::size_t block = 4096, copies = 300;
    std::mt19937 random(718);
    std::string types[2];
    for (std::size_t i = 0; i < block; ++i) types[0].push_back(1+random()%4);
    types[1] = types[0];
    types[1][block-3] = 1+types[1][block-3]%4;
    for (auto& type : types) type.push_back(6);
    std::vector<std::uint8_t> sequence;
    for (std::size_t copy = 0; copy < copies; ++copy)
        sequence.insert(sequence.end(), types[copy%2].begin(), types[copy%2].end());
    auto suffixes = nustar::star2711b::detail::buildSuffixArraySais64(sequence);
    std::erase_if(suffixes, [&](auto position) {
        return position == sequence.size() || sequence[position] == 6;
    });
    for (const std::uint32_t limit : {UINT32_MAX, 1536u}) {
        std::vector<std::pair<std::string, std::uint64_t>> keys;
        for (std::size_t type = 0; type < 2; ++type)
            for (std::size_t offset = 0; offset < block; ++offset)
                keys.emplace_back(types[type].substr(offset, limit), type*block+offset);
        std::sort(keys.begin(), keys.end());
        std::vector<std::uint64_t> expected;
        for (std::size_t first = 0; first < keys.size();) {
            auto last = first+1;
            while (last < keys.size() && keys[last].first == keys[first].first) ++last;
            std::vector<std::uint64_t> positions;
            for (auto key = first; key < last; ++key) {
                const auto type = keys[key].second/block, offset = keys[key].second%block;
                for (auto copy = type; copy < copies; copy += 2)
                    positions.push_back(copy*(block+1)+offset);
            }
            std::sort(positions.begin(), positions.end());
            expected.insert(expected.end(), positions.begin(), positions.end());
            first = last;
        }
        auto actual = suffixes;
        nustar::star2711b::detail::normalizeStarSuffixRanks64(actual, sequence, sequence.size(), 0, limit);
        if (actual != expected) {
            std::cerr << "STAR repeat/finite-prefix rank mismatch: limit=" << limit << '\n';
            return 1;
        }
    }
}
