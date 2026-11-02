#pragma once

// Reserved framework contract; see README for current implementation status.

#include <cstdint>
#include <optional>
#include "../index/native_index.hpp"

namespace nustar::star2711b {

namespace star_compat {
struct JunctionSeedSplit {
    JunctionId junction;
    GenomeOffset donor_start, acceptor_start;
    std::uint64_t donor_length, acceptor_length;
};
Result<std::optional<JunctionSeedSplit>> splitJunctionSeed(
    GenomeOffset synthetic_hit, std::uint64_t seed_length,
    const SyntheticJunctionLayout&, const JunctionDatabaseView&);
} // namespace star_compat

} // namespace nustar::star2711b
