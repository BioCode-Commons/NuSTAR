#pragma once

#include <cstdint>
#include <optional>
#include <limits>
#include <vector>
#include "../common/types.hpp"
#include "../common/stage.hpp"

namespace nustar::star2711b {

// STAR WA_sjA uses UINT64_MAX for an ordinary genomic seed. Keep the
// optional-like API at the module boundary, but store just that one integer.
class WindowJunction {
public:
    static constexpr auto absent = std::numeric_limits<std::uint64_t>::max();
    WindowJunction() = default;
    WindowJunction(std::nullopt_t) {}
    WindowJunction(std::optional<JunctionId> junction)
        : id_(junction.value_or(JunctionId{absent})) {}
    WindowJunction(JunctionId junction) : id_(junction) {}
    bool has_value() const noexcept { return id_.value != absent; }
    explicit operator bool() const noexcept { return has_value(); }
    const JunctionId* operator->() const noexcept { return &id_; }
    std::uint64_t encoded() const noexcept { return id_.value; }
private:
    JunctionId id_{absent};
};

struct WindowSeed {
    ReadOffset read_start; // Coordinates in this window's oriented combined read.
    GenomeOffset genome_start;
    std::uint64_t length, repeat_count;
    bool anchor;
    MateId mate;
    WindowJunction sjdb; // WA_sjA, including STAR's absent sentinel.
};
struct AlignmentWindow {
    WindowId id;
    ChromosomeId chromosome;
    Orientation orientation;
    BinId first_bin, last_bin_inclusive; // STAR WC stores bins, not base positions.
    bool alive;
    std::vector<WindowSeed> seeds; // WA order, including its pruning/tie behavior.
    std::uint64_t minimum_recordable_length;
    std::optional<std::uint64_t> last_anchor;
};
struct WindowBatch {
    ReferenceGeneration generation;
    std::vector<AlignmentWindow> windows;
    StageMarkers markers;
};

} // namespace nustar::star2711b
