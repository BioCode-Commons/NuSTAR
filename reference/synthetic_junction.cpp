#include "synthetic_junction.hpp"

namespace nustar::star2711b {

std::vector<BaseCode> makeSyntheticJunctionSequence(
    std::span<const SjdbEntry> junctions, std::span<const BaseCode> forward,
    std::uint64_t overhang) {
    std::vector<BaseCode> synthetic;
    synthetic.reserve(junctions.size() * (2 * overhang + 1));
    for (const auto& junction : junctions) {
        for (std::uint64_t i = 0; i < overhang; ++i)
            synthetic.push_back(forward[junction.donor_sequence_begin.value + i]);
        for (std::uint64_t i = 0; i < overhang; ++i)
            synthetic.push_back(forward[junction.acceptor_sequence_begin.value + i]);
        synthetic.push_back(BaseCode::GenomePadding);
    }
    return synthetic;
}

} // namespace nustar::star2711b
