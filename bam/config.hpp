#pragma once
#include <cstdint>
namespace nustar::star2711b {
struct BamConfig {
    bool unsorted = false;
    bool coordinate = false;
    int compression = 1;
    std::uint64_t sort_memory = 256ULL << 20;
    unsigned threads = 1;
    bool enabled() const noexcept { return unsorted || coordinate; }
};
} // namespace nustar::star2711b
