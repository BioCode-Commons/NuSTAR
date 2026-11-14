#include "../../seed_search/base_compare_avx2.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/mman.h>
#include <unistd.h>
#endif

using namespace nustar::star2711b;

// Deliberately independent byte oracle: no production comparison helper is used.
template<bool Q, bool G>
std::uint64_t scalar(const BaseCode* query,const BaseCode* genome,
                     std::uint64_t length,bool& greater) {
    for (std::uint64_t i=0;i<length;++i) {
        const auto offset=static_cast<std::ptrdiff_t>(i);
        const auto a=std::uint8_t(query[Q ? offset : -offset]);
        const auto b=std::uint8_t(genome[G ? offset : -offset]);
        if (a!=b) { greater=G ? a>b : !(a>b || b>3);return i; }
    }
    return length;
}

template<bool Q, bool G>
void compare(const BaseCode* query,const BaseCode* genome,std::uint64_t length,
             const BaseCode* begin,const BaseCode* end) {
    for (const bool initial : {false,true}) {
        auto expected_greater=initial,actual_greater=initial;
        const auto expected=scalar<Q,G>(query,genome,length,expected_greater);
        const auto actual=detail::compareBasesAvx2<Q,G>(query,genome,length,begin,end,actual_greater);
        if (actual!=expected || actual_greater!=expected_greater)
            throw std::runtime_error("AVX2 comparison differs from independent byte oracle");
    }
}

template<bool Q, bool G>
void checkDirections() {
    std::array<BaseCode,1024> query{},genome{};
    constexpr std::array lengths{0,1,2,7,8,9,15,16,17,23,24,30,31,32,33,34,
                                47,48,49,63,64,65,95,96,97,127,128,129,150,151,255,256,257};
    // Exercise every unaligned address for AVX2 and every mismatch position,
    // including non-DNA codes and both initial ordering-flag values.
    for (int alignment=0;alignment<32;++alignment) {
        const auto qstart=Q ? 64+alignment : 767+alignment;
        const auto gstart=G ? 96+alignment : 735+alignment;
        for (const auto length : lengths) {
            for (int i=0;i<length;++i) {
                const auto value=BaseCode((i*17+alignment)%4);
                query[qstart+(Q ? i : -i)]=value;
                genome[gstart+(G ? i : -i)]=value;
            }
            const auto check=[&] {
                compare<Q,G>(query.data()+qstart,genome.data()+gstart,length,
                             genome.data(),genome.data()+genome.size());
            };
            check();
            for (int i=0;i<length;++i) {
                auto& a=query[qstart+(Q ? i : -i)];
                auto& b=genome[gstart+(G ? i : -i)];
                const auto saved_a=a,saved_b=b;
                for (unsigned x=0;x<7;++x) for (unsigned y=0;y<7;++y) {
                    a=BaseCode(x);b=BaseCode(y);check();
                }
                a=saved_a;b=saved_b;
            }
            // Boundaries can leave only a short logical genome prefix before
            // scalar guard bytes. Wide reads must stop at the logical boundary.
            for (int available=0;available<=std::min(length,65);++available) {
                const auto* g=genome.data()+gstart;
                const auto* begin=G ? g : g-available+1;
                const auto* end=G ? g+available : g+1;
                compare<Q,G>(query.data()+qstart,g,length,begin,end);
            }
        }
    }
}

#if defined(__unix__) || defined(__APPLE__)
class GuardedPage {
public:
    GuardedPage() : page_size_(static_cast<std::size_t>(::sysconf(_SC_PAGESIZE))) {
        mapping_=::mmap(nullptr,3*page_size_,PROT_NONE,MAP_PRIVATE|MAP_ANON,-1,0);
        if (mapping_==MAP_FAILED) throw std::runtime_error("mmap failed");
        data_=static_cast<BaseCode*>(mapping_)+page_size_;
        if (::mprotect(data_,page_size_,PROT_READ|PROT_WRITE)!=0) {
            ::munmap(mapping_,3*page_size_);
            throw std::runtime_error("mprotect failed");
        }
    }
    ~GuardedPage() { ::munmap(mapping_,3*page_size_); }
    BaseCode* data() const { return data_; }
    std::size_t size() const { return page_size_; }
private:
    std::size_t page_size_;
    void* mapping_;
    BaseCode* data_;
};

template<bool Q, bool G>
void checkGuardPages() {
    GuardedPage query,genome;
    for (int length=1;length<=160;++length) {
        auto* q=Q ? query.data()+query.size()-length : query.data()+length-1;
        auto* g=G ? genome.data()+genome.size()-length : genome.data()+length-1;
        for (int i=0;i<length;++i) {
            q[Q ? i : -i]=BaseCode(i%4);
            g[G ? i : -i]=BaseCode(i%4);
        }
        compare<Q,G>(q,g,length,genome.data(),genome.data()+genome.size());
        // Place a mismatch immediately before the inaccessible guard too.
        g[G ? length-1 : 1-length]=BaseCode::N;
        compare<Q,G>(q,g,length,genome.data(),genome.data()+genome.size());
    }
}
#endif

int main() {
    if (!detail::avx2Available()) {
        std::cout << "AVX2 comparison test skipped: CPU/compiler has no AVX2 kernel\n";
        return 77;
    }
    try {
        checkDirections<true,true>();checkDirections<true,false>();
        checkDirections<false,true>();checkDirections<false,false>();
#if defined(__unix__) || defined(__APPLE__)
        checkGuardPages<true,true>();checkGuardPages<true,false>();
        checkGuardPages<false,true>();checkGuardPages<false,false>();
#endif
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';return 1;
    }
}
