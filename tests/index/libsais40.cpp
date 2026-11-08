#include "../../index/compact_suffix_array.hpp"
#include "../../index/wide_suffix_array.hpp"
#include "../../third_party/libsais40/slot.hpp"
#include "../../third_party/libsais40/native_recursion.hpp"
#include <algorithm>
#include <array>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
using namespace nustar::star2711b;
using namespace nustar::star2711b::detail;
void require(bool b,const char* msg){if(!b)throw std::runtime_error(msg);}
void check(const std::vector<std::uint8_t>& t,int threads,bool naive){
    auto a=buildSuffixArrayLibsais40(t,threads);
    auto b=buildSuffixArraySais64(t,threads);
    require(a.allocatedBytes()==5*(t.size()+1)+7,"workspace is not five bytes per position");
    for(std::size_t i=0;i<b.size();++i)require(a[i]==b[i],"suffix order differs from libsais64");
    if(naive){
        std::vector<std::size_t> p(t.size()+1);std::iota(p.begin(),p.end(),0);
        std::sort(p.begin(),p.end(),[&](auto x,auto y){return std::lexicographical_compare(t.begin()+x,t.end(),t.begin()+y,t.end());});
        require(std::equal(p.begin(),p.end(),b.begin()),"naive suffix order differs");
    }
}
int main(){try{
    using namespace libsais40_impl;
    Slot signed_slots[100];
    for (std::uint64_t bits : {0ULL,1ULL,(1ULL<<31)-1,1ULL<<31,(1ULL<<32)-1,1ULL<<32,(1ULL<<38)-1,1ULL<<38,(1ULL<<39)-1,1ULL<<39,(1ULL<<40)-1}) {
        for (auto& slot:signed_slots) slot=signed40(bits);
        for (const auto& slot:signed_slots) require(unsigned40(std::int64_t(slot))==bits,"signed slot flags or cross-line access failed");
    }
    bool length_rejected=false;
    validateLibsais40Length((1ULL<<38)-2);
    try { validateLibsais40Length((1ULL<<38)-1); }
    catch(const std::length_error&) { length_rejected=true; }
    require(length_rejected,"libsais40 capacity must fail before reading or allocating");
    for (std::size_t n : {64U,4095U,4096U,4097U,8201U}) for (unsigned offset=0;offset<8;++offset) {
        std::vector<Slot> arena(n+offset+8);
        for (auto& slot:arena) slot=1234567;
        auto* source=arena.data()+offset;
        for (std::size_t i=0;i<n;++i) source[i]=(i*104729ULL)%INT32_MAX;
        auto* narrowed=narrowInPlace(source,n);
        require(reinterpret_cast<std::uintptr_t>(narrowed)%alignof(std::int32_t)==0,"recursive input is unaligned");
        for (std::size_t i=0;i<n;++i) require(narrowed[i]==std::int64_t((i*104729ULL)%INT32_MAX),"narrowing overwrote unread slots");
        expandInPlace(source,narrowed,n);
        for (std::size_t i=0;i<n;++i) require(std::int64_t(source[i])==std::int64_t((i*104729ULL)%INT32_MAX),"expansion overwrote unread slots");
        for (unsigned i=0;i<offset;++i) require(std::int64_t(arena[i])==1234567,"recursive prefix guard overwritten");
        for (std::size_t i=offset+n;i<arena.size();++i) require(std::int64_t(arena[i])==1234567,"recursive suffix guard overwritten");
    }
    CompactSuffixArray s(100);
    const std::array<std::uint64_t,8> values{0,1,(1ULL<<31)-1,1ULL<<31,(1ULL<<32)-1,1ULL<<32,Position40::flag-1,Position40::mask};
    for(std::size_t i=0;i<s.size();++i)s[i]=values[i%values.size()];
    for(std::size_t i=0;i<s.size();++i)require(s[i]==values[i%values.size()],"position or state bit truncated");
    bool rejected=false;try{CompactSuffixArray bad(Position40::flag);}catch(const std::length_error&){rejected=true;}
    require(rejected,"capacity must fail before allocation");
    for(unsigned n=0;n<=7;++n){unsigned combinations=1;for(unsigned i=0;i<n;++i)combinations*=3;
        for(unsigned code=0;code<combinations;++code){auto x=code;std::vector<std::uint8_t> t(n);for(auto& c:t){c=1+x%3;x/=3;}check(t,1,true);}}
    std::mt19937 rng(271140);
    for(int n=0;n<1024;n+=13)for(int alphabet=1;alphabet<=6;++alphabet){
        std::vector<std::uint8_t> t(n);for(auto& c:t)c=1+rng()%alphabet;check(t,1,true);}
    // Exercise parallel induction and recursive integer alphabets.
    for(int kind=0;kind<4;++kind){std::vector<std::uint8_t> t((1<<20)+37);
        for(std::size_t i=0;i<t.size();++i)t[i]=kind==0?1:kind==1?1+i%2:kind==2?1+rng()%6:1+(i%4093)%5;
        for(int threads:{1,2,4,32})check(t,threads,false);}
    // Cover lengths around powers of two and uneven OpenMP partitions.
    for(std::size_t n:{1U<<20,(1U<<20)+1,(2U<<20)-1,2U<<20}) {
        std::vector<std::uint8_t> t(n);
        for(auto& c:t)c=1+rng()%6;
        for(int threads:{3,32})check(t,threads,false);
    }
    std::cout<<"libsais40 storage, exhaustive, random, repetitive and parallel comparisons passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
