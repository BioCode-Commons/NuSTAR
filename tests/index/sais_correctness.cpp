#include "../../index/sais.hpp"
#include "../../index/wide_suffix_array.hpp"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

using nustar::star2711b::buildSuffixArraySais;

int main() {
    std::mt19937 random(2711);
    for (int length = 0; length <= 128; ++length) {
        for (int trial = 0; trial < 80; ++trial) {
            std::vector<std::uint8_t> sequence(length);
            const auto alphabet = 1 + trial % 6;
            for (auto& symbol : sequence) symbol = 1 + random() % alphabet;
            auto actual = buildSuffixArraySais(sequence);
            std::vector<std::uint32_t> expected(sequence.size()+1);
            std::iota(expected.begin(), expected.end(), 0);
            std::sort(expected.begin(), expected.end(), [&](auto a, auto b) {
                return std::lexicographical_compare(sequence.begin()+a, sequence.end(),
                                                    sequence.begin()+b, sequence.end());
            });
            if (actual != expected) {
                std::cerr << "SAIS mismatch at length " << length
                          << ", trial " << trial << '\n';
                return 1;
            }
            const auto wide = nustar::star2711b::detail::buildSuffixArraySais64(sequence);
            if (wide.size() != expected.size() || !std::equal(wide.begin(), wide.end(), expected.begin())) {
                std::cerr << "64-bit SAIS mismatch at length " << length << ", trial " << trial << '\n';
                return 1;
            }
        }
    }
    return 0;
}
