#include "../../seed_search/seed_compare.hpp"

#include <array>
#include <iostream>
#include <stdexcept>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/mman.h>
#include <unistd.h>
#endif

using namespace nustar::star2711b;

template<bool Q, bool G>
void checkComparisons(const BaseCode* query, const BaseCode* genome,
    std::uint64_t length, const BaseCode* begin, const BaseCode* end,
    std::uint64_t expected, bool expected_greater, bool initial) {
    bool greater = initial;
    const auto actual = detail::compareBasesNeon<Q,G>(query, genome, length,
        begin, end, greater);
    if (actual != expected || greater != expected_greater)
        throw std::runtime_error("NEON comparison differs from STAR byte ordering");
    for (const auto use_avx2 : {false, true}) {
        if (use_avx2 && !detail::avx2Available()) continue;
        greater = initial;
        const auto seed = detail::compareSeedBases<Q,G>(query, genome, length,
            begin, end, greater, use_avx2);
        if (seed != expected || greater != expected_greater)
            throw std::runtime_error("seed comparison differs from STAR byte ordering");
    }
}

template<bool Q, bool G>
void checkDirections() {
    std::array<BaseCode, 768> query{}, genome{};
    constexpr std::array lengths{0,1,2,3,7,8,9,14,15,16,17,18,23,24,25,
                                30,31,32,33,47,48,49,63,64,65,95,96,97,
                                127,128,129,150,151,160,161,255,256,257};
    // Logical bounds deliberately differ from allocated bounds, exercising
    // scalar guard comparisons and the transition into/out of vector loads.
    for (const auto start : {0,1,14,15,16,19,20,21,27,28,31,32,33,100,
                            251,252,259,260,264,265,279,280,281,300,751,752,766,767}) {
        for (const auto requested : lengths) {
            const auto length = std::min(requested, G ? 768-start : start+1);
            const auto qstart = Q ? 37+requested%16 : 600+requested%16;
            for (int i = 0; i < length; ++i) {
                const auto base = BaseCode((i*17+start)%4);
                query[qstart+(Q ? i : -i)] = base;
                genome[start+(G ? i : -i)] = base;
            }
            const auto check = [&] {
                for (const auto initial : {false, true}) {
                    bool expected_greater = initial;
                    std::uint64_t expected = 0;
                    // Independent byte oracle: do not compare two SIMD paths.
                    for (; expected < std::uint64_t(length); ++expected) {
                        const auto i = static_cast<int>(expected);
                        const auto a = std::uint8_t(query[qstart+(Q ? i : -i)]);
                        const auto b = std::uint8_t(genome[start+(G ? i : -i)]);
                        if (a != b) {
                            expected_greater = G ? a > b : !(a > b || b > 3);
                            break;
                        }
                    }
                    checkComparisons<Q,G>(query.data()+qstart,
                        genome.data()+start, length, genome.data()+20, genome.data()+280,
                        expected, expected_greater, initial);
                }
            };
            check();
            for (int i = 0; i < length; ++i) {
                auto& q = query[qstart+(Q ? i : -i)];
                auto& g = genome[start+(G ? i : -i)];
                const auto saved_q = q, saved_g = g;
                for (unsigned a = 0; a < 6; ++a) for (unsigned b = 0; b < 6; ++b) {
                    q = BaseCode(a); g = BaseCode(b); check();
                }
                q = saved_q; g = saved_g;
            }
        }
    }
}

#if defined(__unix__) || defined(__APPLE__)
class GuardedPage {
public:
    GuardedPage() : size_(static_cast<std::size_t>(sysconf(_SC_PAGESIZE))) {
        mapping_ = mmap(nullptr, size_*3, PROT_NONE, MAP_PRIVATE|MAP_ANON, -1, 0);
        if (mapping_ == MAP_FAILED) throw std::runtime_error("mmap failed");
        data_ = reinterpret_cast<BaseCode*>(mapping_)+size_;
        if (mprotect(data_, size_, PROT_READ|PROT_WRITE)) {
            munmap(mapping_, size_*3);
            throw std::runtime_error("mprotect failed");
        }
    }
    ~GuardedPage() { munmap(mapping_, size_*3); }
    BaseCode* begin() const { return data_; }
    BaseCode* end() const { return data_+size_; }
private:
    std::size_t size_;
    void* mapping_;
    BaseCode* data_;
};

template<bool Q, bool G>
void checkGuardPages() {
    GuardedPage query, genome;
    // Match all the way up to an inaccessible page in either direction. Any
    // speculative vector/word read outside the requested span faults here.
    for (int length = 1; length <= 257; ++length) {
        auto* q = Q ? query.end()-length : query.begin()+length-1;
        auto* g = G ? genome.end()-length : genome.begin()+length-1;
        for (int i = 0; i < length; ++i) {
            q[Q ? i : -i] = BaseCode(i%4);
            g[G ? i : -i] = BaseCode(i%4);
        }
        for (bool initial : {false, true}) {
            checkComparisons<Q,G>(q, g, length, genome.begin(), genome.end(),
                length, initial, initial);
        }
        // Also require the first unequal byte to win over later mismatches.
        const auto saved_last = g[G ? length-1 : 1-length];
        g[G ? length-1 : 1-length] = BaseCode::N;
        for (int first = 0; first < length; ++first) {
            const auto saved = g[G ? first : -first];
            g[G ? first : -first] = BaseCode::N;
            for (bool initial : {false, true})
                checkComparisons<Q,G>(q, g, length, genome.begin(), genome.end(),
                    first, false, initial);
            g[G ? first : -first] = saved;
        }
        g[G ? length-1 : 1-length] = saved_last;
    }
}
#endif

int main() {
    try {
        checkDirections<true,true>(); checkDirections<true,false>();
        checkDirections<false,true>(); checkDirections<false,false>();
#if defined(__unix__) || defined(__APPLE__)
        checkGuardPages<true,true>(); checkGuardPages<true,false>();
        checkGuardPages<false,true>(); checkGuardPages<false,false>();
#endif
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
