#include "compact_suffix_array.hpp"
#include "wide_suffix_array.hpp"
#include <iostream>
#include <new>
#include <stdexcept>
#include <string>
#include "../third_party/libsais40/libsais40.inc"

namespace nustar::star2711b::detail {
void validateLibsais40Length(std::size_t length) {
    if (length>=(std::uint64_t{1}<<38)-1)
        throw std::length_error("libsais40 requires physical length < 2^38-1");
}
CompactSuffixArray buildSuffixArrayLibsais40(std::span<const std::uint8_t> text,
                                           int runThreadN) {
    // Leave both the sign marker and the suffix-group marker available.
    validateLibsais40Length(text.size());
    CompactSuffixArray result(text.size()+1);
    result[0]=text.size();
    if (text.empty()) return result;
    if (text.size()==1) { result[1]=0; return result; }
    const auto threads=suffixConstructionThreads(text.size(),runThreadN);
    using libsais40_impl::Slot;
    auto* storage=result.begin()+1;
    auto* signed_slots=reinterpret_cast<Slot*>(storage);
    // Start the temporary signed-slot lifetimes in the same raw storage.
    // Trivial default initialization preserves bytes and emits no clearing pass.
    for (std::size_t i=0;i<text.size();++i)
        ::new (static_cast<void*>(signed_slots+i)) Slot;
    const auto status=libsais40_impl::libsais40_main(text.data(),signed_slots,
        static_cast<std::int64_t>(text.size()),0,0,nullptr,0,nullptr,threads);
    // All surviving entries are nonnegative positions. Restore the unsigned
    // public slot type without moving or expanding the suffix array.
    for (std::size_t i=0;i<text.size();++i)
        ::new (static_cast<void*>(storage+i)) Position40;
    if (status!=0)
        throw std::runtime_error("libsais40 suffix sorting failed: "+std::to_string(status));
    std::cerr << "Suffix construction: libsais40-v1 (libsais 2.12.1), positions=40bit, threads="
              << threads << ", physical_bases=" << text.size()
              << ", sa_workspace_bytes=" << result.allocatedBytes() << '\n';
    return result;
}
}
