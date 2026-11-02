#include "sais.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace nustar::star2711b {
namespace {

// Induced sorting, including the unique zero sentinel. The reduced LMS string is
// recursively sorted by the same procedure. This does not rely on NuSTAR's packed
// array allocation or its in-place sentinel arithmetic.
template<class Symbol>
std::vector<std::uint32_t> sais(std::span<const Symbol> symbols, std::int32_t upper) {
    constexpr auto vacant = std::numeric_limits<std::uint32_t>::max();
    const auto n = static_cast<std::int32_t>(symbols.size());
    if (n == 0) return {};
    if (n == 1) return {0};
    if (n == 2) return symbols[0] < symbols[1] ?
        std::vector<std::uint32_t>{0,1} : std::vector<std::uint32_t>{1,0};
    std::vector<std::uint32_t> sa(n, vacant);
    std::vector<std::uint8_t> is_s(n, 0);
    is_s[n-1] = 1;
    for (std::int32_t i = n-2; i >= 0; --i)
        is_s[i] = symbols[i] == symbols[i+1] ? is_s[i+1] : symbols[i] < symbols[i+1];
    std::vector<std::int32_t> sum_l(upper+2), sum_s(upper+2);
    for (std::int32_t i = 0; i < n; ++i) {
        if (is_s[i]) ++sum_l[symbols[i]+1];
        else ++sum_s[symbols[i]];
    }
    for (std::int32_t i = 0; i <= upper; ++i) {
        sum_s[i] += sum_l[i];
        sum_l[i+1] += sum_s[i];
    }
    auto induce = [&](std::span<const std::int32_t> lms) {
        std::fill(sa.begin(), sa.end(), vacant);
        auto buf = sum_s;
        for (const auto position : lms)
            if (position < n) sa[buf[symbols[position]]++] = position;
        buf = sum_l;
        sa[buf[symbols[n-1]]++] = n-1;
        for (std::int32_t i = 0; i < n; ++i) {
            const auto position = sa[i];
            if (position != vacant && position >= 1 && !is_s[position-1])
                sa[buf[symbols[position-1]]++] = position-1;
        }
        buf = sum_l;
        for (std::int32_t i = n-1; i >= 0; --i) {
            const auto position = sa[i];
            if (position != vacant && position >= 1 && is_s[position-1])
                sa[--buf[symbols[position-1]+1]] = position-1;
        }
    };
    std::vector<std::int32_t> lms_map(n+1, -1), lms;
    for (std::int32_t i = 1; i < n; ++i)
        if (!is_s[i-1] && is_s[i]) {
            lms_map[i] = static_cast<std::int32_t>(lms.size());
            lms.push_back(i);
        }
    induce(lms);
    if (!lms.empty()) {
        std::vector<std::int32_t> ordered;
        ordered.reserve(lms.size());
        for (const auto position : sa)
            if (position != vacant && lms_map[position] >= 0) ordered.push_back(position);
        std::vector<std::int32_t> reduced(lms.size());
        std::int32_t upper_reduced = 0;
        reduced[lms_map[ordered[0]]] = 0;
        for (std::size_t i = 1; i < ordered.size(); ++i) {
            const auto left = ordered[i-1], right = ordered[i];
            const auto left_end = lms_map[left]+1 < static_cast<std::int32_t>(lms.size()) ?
                lms[lms_map[left]+1] : n;
            const auto right_end = lms_map[right]+1 < static_cast<std::int32_t>(lms.size()) ?
                lms[lms_map[right]+1] : n;
            bool same = left_end-left == right_end-right;
            if (same) {
                for (std::int32_t offset = 0; offset <= left_end-left; ++offset) {
                    const auto a = left+offset, b = right+offset;
                    if (a >= n || b >= n) {
                        same = a >= n && b >= n;
                        break;
                    }
                    if (symbols[a] != symbols[b] || is_s[a] != is_s[b]) {
                        same = false;
                        break;
                    }
                }
            }
            if (!same) ++upper_reduced;
            reduced[lms_map[right]] = upper_reduced;
        }
        const auto reduced_sa = sais(std::span<const std::int32_t>(reduced), upper_reduced);
        for (std::size_t i = 0; i < ordered.size(); ++i)
            ordered[i] = lms[reduced_sa[i]];
        induce(ordered);
    }
    return sa;
}

} // namespace

std::vector<std::uint32_t> buildSuffixArraySais(std::span<const std::uint8_t> symbols) {
    if (symbols.size() >= static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
        throw std::length_error("SAIS input exceeds 32-bit index capacity");
    std::vector<std::uint8_t> input;
    input.reserve(symbols.size()+1);
    std::int32_t upper = 0;
    for (const auto symbol : symbols) {
        if (symbol == 0) throw std::invalid_argument("SAIS input contains sentinel");
        input.push_back(symbol);
        upper = std::max(upper, static_cast<std::int32_t>(symbol));
    }
    input.push_back(0);
    return sais(std::span<const std::uint8_t>(input), upper);
}

} // namespace nustar::star2711b
