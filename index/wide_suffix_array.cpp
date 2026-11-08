#include "wide_suffix_array.hpp"
#include "../third_party/libsais/libsais64.h"

#include <algorithm>
#include <charconv>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <thread>

namespace nustar::star2711b::detail {

int suffixConstructionThreads(std::size_t length, int runThreadN) {
#if defined(LIBSAIS_OPENMP)
    int threads = static_cast<int>(std::min(16u, std::max(1u, std::thread::hardware_concurrency())));
    if (runThreadN > 0) threads = runThreadN;
    else if (const auto* configured = std::getenv("NUSTAR_SA_THREADS")) {
        const std::string_view text(configured);
        const auto parsed = std::from_chars(text.data(), text.data()+text.size(), threads);
        if (parsed.ec != std::errc{} || parsed.ptr != text.data()+text.size() || threads < 1 || threads > 256)
            throw std::invalid_argument("NUSTAR_SA_THREADS must be between 1 and 256");
    }
    return length < (1ULL << 20) ? 1 : threads;
#else
    (void)length;
    (void)runThreadN;
    return 1;
#endif
}

std::vector<std::uint64_t> buildSuffixArraySais64(std::span<const std::uint8_t> symbols, int runThreadN) {
    if (symbols.size() >= static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
        throw std::length_error("suffix array input exceeds signed 64-bit capacity");
    std::vector<std::uint64_t> suffixes(symbols.size()+1);
    // libsais has an implicit end marker. Prepending its empty suffix yields
    // exactly the existing builder contract without copying a 6+ GB input.
    suffixes[0] = symbols.size();
    if (symbols.empty()) return suffixes;
    const auto threads = suffixConstructionThreads(symbols.size(), runThreadN);
    std::int64_t status;
#if defined(LIBSAIS_OPENMP)
    status = libsais64_omp(symbols.data(), reinterpret_cast<std::int64_t*>(suffixes.data()+1),
                          static_cast<std::int64_t>(symbols.size()), 0, nullptr, threads);
#else
    status = libsais64(symbols.data(), reinterpret_cast<std::int64_t*>(suffixes.data()+1),
                      static_cast<std::int64_t>(symbols.size()), 0, nullptr);
#endif
    if (status != 0) throw std::runtime_error("libsais64 suffix sorting failed: " + std::to_string(status));
    std::cerr << "Suffix construction: libsais64 2.12.1, positions=64bit, threads=" << threads
              << ", physical_bases=" << symbols.size() << '\n';
    return suffixes;
}

} // namespace nustar::star2711b::detail
