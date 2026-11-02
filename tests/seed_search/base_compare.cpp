#include "../../seed_search/base_compare.hpp"

#include <array>
#include <iostream>
#include <stdexcept>

using namespace nustar::star2711b;

template<bool Q, bool G>
void checkDirections() {
    std::array<BaseCode,512> query{},genome{};
    constexpr std::array lengths{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,
                                31,32,33,63,64,65,127,128,129,150,151,160,161};
    for (const auto start : {0,1,7,8,19,20,21,27,28,100,251,252,259,260,279,280,281,300,504,505,510,511}) {
        for (const auto requested : lengths) {
            const auto length=std::min(requested,G ? 512-start : start+1);
            const auto qstart=Q ? 37+requested%8 : 400+requested%8;
            for (int i=0;i<length;++i) {
                const auto base=BaseCode((i*17+start)%4);
                query[qstart+(Q ? i : -i)]=base;
                genome[start+(G ? i : -i)]=base;
            }
            const auto check=[&] {
                for (const auto initial : {false,true}) {
                    bool expected_greater=initial,actual_greater=initial;
                    std::uint64_t expected=0;
                    // Independent STAR-style byte loop, including its reverse
                    // genomic N ordering and unchanged flag on full equality.
                    for (;expected<std::uint64_t(length);++expected) {
                        const int i=static_cast<int>(expected);
                        const auto a=std::uint8_t(query[qstart+(Q ? i : -i)]);
                        const auto b=std::uint8_t(genome[start+(G ? i : -i)]);
                        if (a!=b) { expected_greater=G ? a>b : !(a>b || b>3); break; }
                    }
                    const auto actual=detail::compareBases<Q,G>(query.data()+qstart,
                        genome.data()+start,length,genome.data()+20,genome.data()+280,
                        actual_greater);
                    if (actual!=expected || actual_greater!=expected_greater)
                        throw std::runtime_error("word comparison differs from STAR byte ordering");
                }
            };
            check();
            for (int i=0;i<length;++i) {
                auto& q=query[qstart+(Q ? i : -i)];
                auto& g=genome[start+(G ? i : -i)];
                const auto saved_q=q,saved_g=g;
                for (unsigned a=0;a<6;++a) for (unsigned b=0;b<6;++b) {
                    q=BaseCode(a);g=BaseCode(b);check();
                }
                q=saved_q;g=saved_g;
            }
        }
    }
}

int main() {
    try {
        checkDirections<true,true>();checkDirections<true,false>();
        checkDirections<false,true>();checkDirections<false,false>();
    } catch (const std::exception& e) { std::cerr << e.what() << '\n';return 1; }
}
