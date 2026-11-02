#include "star_rank.hpp"
#include "wide_suffix_array.hpp"

#include <algorithm>
#include <stdexcept>
#include <unordered_map>

namespace nustar::star2711b {
namespace {

struct CommonInterval {
    std::uint64_t begin, end;
    bool terminated;
};

// An exact, bounded cache of already compared intervals, keyed by displacement.
// Unlike a rolling hash, this cannot treat different strings as equal. Long
// duplicated chromosome segments would otherwise be reread at every SA rank.
using ComparisonCache = std::unordered_map<std::uint64_t, CommonInterval>;

bool sameStarSuffix(std::span<const std::uint8_t> sequence,
                    std::uint64_t left, std::uint64_t right,
                    std::uint32_t suffix_length_max, ComparisonCache* cache = nullptr) {
    if (left > right) std::swap(left, right);
    const auto displacement = right-left;
    if (cache) {
        const auto found = cache->find(displacement);
        if (found != cache->end() && left >= found->second.begin && left < found->second.end)
            return found->second.terminated || suffix_length_max <= found->second.end-left;
    }
    std::uint64_t i = 0;
    bool terminated = false;
    for (; i < suffix_length_max; ++i) {
        if (left+i >= sequence.size() || right+i >= sequence.size()) break;
        const auto a = sequence[left+i], b = sequence[right+i];
        if (a != b) break;
        if (a == 6) { terminated = true; break; }
    }
    if (cache && i >= 1024 && i < suffix_length_max) {
        if (cache->size() >= 8192) cache->clear();
        (*cache)[displacement] = CommonInterval{left, left+i, terminated};
    }
    return terminated || i == suffix_length_max;
}

std::uint64_t starCoordinate(std::uint64_t position, std::uint64_t real,
                             std::uint64_t synthetic) {
    const auto combined = real+synthetic;
    const auto forward_synthetic = 2*real+40;
    const auto reverse_synthetic = forward_synthetic+synthetic+20;
    if (position < real) return position;
    if (position >= real+20 && position < 2*real+20)
        return combined+synthetic+position-(real+20);
    if (position >= forward_synthetic && position < forward_synthetic+synthetic)
        return real+position-forward_synthetic;
    if (position >= reverse_synthetic && position < reverse_synthetic+synthetic)
        return combined+position-reverse_synthetic;
    throw std::logic_error("SA contains a separator position");
}

template<class Position>
void normalize(std::span<Position> suffixes,
    std::span<const std::uint8_t> sequence, std::uint64_t real_length,
    std::uint64_t synthetic_length, std::uint32_t suffix_length_max, int runThreadN = 0) {
    if (suffix_length_max == 0) throw std::invalid_argument("empty STAR suffix comparison");
    std::vector<std::uint8_t> equal;
    if (suffixes.size() >= (1ULL << 20)) {
        // STAR equality ends at padding. Equal terminated suffixes necessarily
        // have equal distances to padding, hence equal distances modulo 256.
        // This inexpensive filter avoids rereading full LCPs for billions of
        // unequal neighbours, while keeping exact comparisons for possible ties.
        std::vector<std::uint8_t> distance(sequence.size());
        std::uint64_t run = 0, longest = 0;
        for (std::size_t i = sequence.size(); i-- > 0;) {
            run = sequence[i] == 6 ? 0 : run+1;
            longest = std::max(longest, run);
            distance[i] = static_cast<std::uint8_t>(run);
        }
        const bool filter = longest < suffix_length_max;
        equal.resize(suffixes.size(), 0);
        const auto threads = detail::suffixConstructionThreads(sequence.size(), runThreadN);
#if defined(_OPENMP)
        #pragma omp parallel num_threads(threads)
#endif
        {
            ComparisonCache cache;
#if defined(_OPENMP)
            #pragma omp for schedule(static)
#endif
            for (std::int64_t i = 1; i < static_cast<std::int64_t>(suffixes.size()); ++i) {
                const auto left = suffixes[i-1], right = suffixes[i];
                equal[i] = (!filter || distance[left] == distance[right]) &&
                    sameStarSuffix(sequence, left, right, suffix_length_max, &cache);
            }
        }
        (void)threads;
    }
    for (std::size_t first = 0; first < suffixes.size();) {
        auto last = first+1;
        while (last < suffixes.size() &&
               (equal.empty() ? sameStarSuffix(sequence, suffixes[last-1], suffixes[last], suffix_length_max)
                              : equal[last] != 0))
            ++last;
        if (last-first > 1)
            std::sort(suffixes.begin()+first, suffixes.begin()+last,
                [&](auto a, auto b) {
                    return starCoordinate(a, real_length, synthetic_length) <
                           starCoordinate(b, real_length, synthetic_length);
                });
        first = last;
    }
}

} // namespace

void normalizeStarSuffixRanks(std::span<std::uint32_t> suffixes,
    std::span<const std::uint8_t> sequence, std::uint64_t real_length,
    std::uint64_t synthetic_length, std::uint32_t suffix_length_max) {
    normalize(suffixes, sequence, real_length, synthetic_length, suffix_length_max);
}

void detail::normalizeStarSuffixRanks64(std::span<std::uint64_t> suffixes,
    std::span<const std::uint8_t> sequence, std::uint64_t real_length,
    std::uint64_t synthetic_length, std::uint32_t suffix_length_max, int runThreadN) {
    normalize(suffixes, sequence, real_length, synthetic_length, suffix_length_max, runThreadN);
}

} // namespace nustar::star2711b
