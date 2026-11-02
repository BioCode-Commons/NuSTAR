#include "physical_sequence.hpp"

#include <stdexcept>
#include <limits>

namespace nustar::star2711b {

std::vector<std::uint8_t> makeStarPhysicalSequence(
    const GenomeReferenceData& reference, std::uint8_t shift) {
    if (shift > 1) throw std::invalid_argument("physical sequence shift must be 0 or 1");
    const auto& forward = reference.storage().forward;
    const auto synthetic = reference.synthetic();
    constexpr auto maximum = (std::numeric_limits<std::size_t>::max()-80)/2;
    if (forward.size() > maximum || synthetic.size() > maximum-forward.size())
        throw std::length_error("physical sequence exceeds addressable memory");
    std::vector<std::uint8_t> sequence;
    sequence.reserve(2*(forward.size()+synthetic.size())+80);
    auto append = [&](std::span<const BaseCode> bases) {
        for (auto base : bases)
            sequence.push_back(static_cast<std::uint8_t>(base)+shift);
    };
    auto appendReverse = [&](std::span<const BaseCode> bases) {
        for (auto it = bases.rbegin(); it != bases.rend(); ++it) {
            const auto value = static_cast<std::uint8_t>(*it);
            sequence.push_back((value < 4 ? (value^3) : value)+shift);
        }
    };
    auto gap = [&] { sequence.insert(sequence.end(), 20, 5+shift); };
    append(forward); gap();
    appendReverse(forward); gap();
    append(synthetic); gap();
    appendReverse(synthetic); gap();
    return sequence;
}

} // namespace nustar::star2711b
