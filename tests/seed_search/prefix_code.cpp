#include "../../seed_search/prefix_code.hpp"

#include <array>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/mman.h>
#include <unistd.h>
#endif

using namespace nustar::star2711b;

std::uint64_t oracle(const BaseCode* query,unsigned length) {
    std::uint64_t result=0;
    for (unsigned i=0;i<length;++i) result=result*4|std::uint8_t(query[i]);
    return result;
}

void check(const BaseCode* query,unsigned length) {
    if (detail::prefixCode(query,length)!=oracle(query,length))
        throw std::runtime_error("prefix encoder differs from scalar shift/OR semantics");
}

int main() {
    try {
        std::array<BaseCode,48> query{};
        // Exhaustively enumerate every canonical sequence through length ten.
        for (unsigned length=0;length<=10;++length) {
            const auto total=std::uint64_t{1}<<(2*length);
            for (std::uint64_t code=0;code<total;++code) {
                auto value=code;
                for (unsigned i=length;i>0;--i) {
                    query[i-1]=BaseCode(value&3);value>>=2;
                }
                check(query.data(),length);
                if (detail::prefixCode(query.data(),length)!=code)
                    throw std::runtime_error("canonical prefix has incorrect numeric value");
            }
        }
        // Every unaligned start, long prefixes, and arbitrary byte codes.
        std::mt19937_64 rng(2711);
        for (unsigned trial=0;trial<200000;++trial) {
            const auto length=unsigned(rng()%16);
            const auto start=unsigned(rng()%16);
            for (unsigned i=0;i<length;++i)
                query[start+i]=BaseCode(rng()%(trial%3==0 ? 256 : trial%3==1 ? 7 : 4));
            check(query.data()+start,length);
        }
        // Place each noncanonical byte at every position, including the tail.
        for (unsigned length=1;length<=15;++length) {
            for (unsigned i=0;i<length;++i) query[i]=BaseCode(i%4);
            for (unsigned i=0;i<length;++i) {
                const auto saved=query[i];
                for (unsigned code=4;code<256;++code) {
                    query[i]=BaseCode(code);check(query.data(),length);
                }
                query[i]=saved;
            }
        }
#if defined(__unix__) || defined(__APPLE__)
        const auto page=static_cast<std::size_t>(::sysconf(_SC_PAGESIZE));
        auto* memory=::mmap(nullptr,page*3,PROT_NONE,MAP_PRIVATE|MAP_ANON,-1,0);
        if (memory==MAP_FAILED) throw std::runtime_error("mmap failed");
        auto* begin=static_cast<BaseCode*>(memory)+page;
        if (::mprotect(begin,page,PROT_READ|PROT_WRITE)!=0) {
            ::munmap(memory,page*3);throw std::runtime_error("mprotect failed");
        }
        for (unsigned length=0;length<=15;++length) {
            auto* query=begin+page-length;
            for (unsigned i=0;i<length;++i) query[i]=BaseCode(i%4);
            check(query,length); // Empty span points directly at a guard page.
        }
        ::munmap(memory,page*3);
#endif
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';return 1;
    }
}
