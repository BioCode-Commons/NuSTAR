#pragma once
#include "slot.hpp"
#include "../libsais/libsais.h"
#include <algorithm>
#include <array>
#include <climits>
#include <new>

namespace nustar::star2711b::detail::libsais40_impl {
inline std::int32_t* narrowInPlace(Slot* source, std::size_t n) {
    const auto address=(reinterpret_cast<std::uintptr_t>(source)+3)&~std::uintptr_t{3};
    auto* target=reinterpret_cast<std::int32_t*>(address);
    std::array<std::int32_t,4096> buffer;
    for (std::size_t first=0;first<n;) {
        const auto count=std::min(buffer.size(),n-first);
        for (std::size_t i=0;i<count;++i) buffer[i]=std::int64_t(source[first+i]);
        for (std::size_t i=0;i<count;++i) ::new (static_cast<void*>(target+first+i)) std::int32_t(buffer[i]);
        first+=count;
    }
    return target;
}
inline void expandInPlace(Slot* target,const std::int32_t* source,std::size_t n) {
    std::array<std::int32_t,4100> buffer;
    while (n) {
        auto count=std::min(std::size_t{4096},n);
        // A shifted aligned source can overlap the first three destination
        // slots. Consume that final short prefix in the same buffered block.
        if (n-count<4) count=n;
        const auto first=n-count;
        for (std::size_t i=0;i<count;++i) buffer[i]=source[first+i];
        for (std::size_t i=0;i<count;++i) {
            auto* slot=::new (static_cast<void*>(target+first+i)) Slot;
            *slot=buffer[i];
        }
        n=first;
    }
}
// Preserve the parent's free-space boundary, including any compacted unique LMS
// entries beyond it. Alignment costs at most three bytes. No new SA is allocated.
inline std::int64_t nativeRecursion(Slot* text,Slot* sa,std::int64_t n,
                                    std::int64_t k,std::int64_t fs,int threads) {
    auto* input=narrowInPlace(text,n);
    const auto address=(reinterpret_cast<std::uintptr_t>(sa)+3)&~std::uintptr_t{3};
    auto* output=reinterpret_cast<std::int32_t*>(address);
    const auto alignment=address-reinterpret_cast<std::uintptr_t>(sa);
    const auto capacity=((std::uint64_t(n)+fs)*5-alignment)/4;
    const auto free_slots=std::min<std::uint64_t>(capacity-n,INT32_MAX-n);
#if defined(LIBSAIS_OPENMP)
    const auto status=libsais_int_omp(input,output,n,k,free_slots,threads);
#else
    (void)threads;
    const auto status=libsais_int(input,output,n,k,free_slots);
#endif
    if (status==0) {
        expandInPlace(text,input,n);
        expandInPlace(sa,output,n);
    }
    return status;
}
}
