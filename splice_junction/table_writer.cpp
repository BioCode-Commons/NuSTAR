#include "table_writer.hpp"
#include <algorithm>
#include <sstream>

namespace nustar::star2711b {

namespace {
std::vector<std::byte> bytes(std::string_view text) {
    std::vector<std::byte> result;
    result.reserve(text.size());
    for (unsigned char character : text)
        result.push_back(static_cast<std::byte>(character));
    return result;
}

} // namespace

Result<std::vector<std::byte>> encodeSjTable(const FinalJunctions& junctions,
                                             const ReferenceMetadata& reference) {
    std::ostringstream output;
    for (const auto& junction : junctions.accepted) {
        auto after = std::upper_bound(reference.chromosomes.begin(),
            reference.chromosomes.end(), junction.key.start.value,
            [](std::uint64_t value, const Chromosome& chromosome) {
                return value < chromosome.start.value;
            });
        const auto& chromosome = *std::prev(after);
        output << chromosome.name << '\t'
               << junction.key.start.value+1-chromosome.start.value << '\t'
               << junction.key.start.value+junction.key.gap-
                      chromosome.start.value << '\t'
               << unsigned(junction.strand) << '\t'
               << unsigned(junction.motif) << '\t'
               << unsigned(junction.annotation) << '\t'
               << junction.unique_count << '\t' << junction.multiple_count << '\t'
               << junction.max_left_overhang << '\n';
    }
    return bytes(output.str());
}

} // namespace nustar::star2711b
