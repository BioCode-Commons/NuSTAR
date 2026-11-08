#include "../../index/packed_suffix_array.hpp"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <ranges>
#include <stdexcept>
#include <vector>

using nustar::star2711b::PackedSuffixArray;

int main() {
    try {
        for (const std::uint8_t width :
             {1, 2, 7, 8, 9, 15, 16, 25, 31, 32, 33, 34, 40, 57, 63, 64}) {
            const std::uint64_t physical = width == 64 ?
                std::numeric_limits<std::uint64_t>::max() :
                std::uint64_t{1} << width;
            std::vector<std::uint64_t> values(257);
            for (std::size_t i = 0; i < values.size(); ++i)
                values[i] = i%3 == 0 ? physical-1 : (i*104729)%physical;
            const PackedSuffixArray packed(values, physical);
            std::vector<std::uint8_t> expected((values.size()*width+7)/8, 0);
            for (std::size_t rank = 0; rank < values.size(); ++rank)
                for (unsigned bit = 0; bit < width; ++bit)
                    if ((values[rank] >> bit) & 1)
                        expected[(rank*width+bit)/8] |=
                            std::uint8_t{1} << ((rank*width+bit)%8);
            if (packed.bitWidth() != width ||
                packed.packedBytes().size() != (values.size()*width+7)/8 ||
                !std::ranges::equal(packed.packedBytes(), expected) ||
                !std::ranges::equal(packed, values))
                throw std::runtime_error("packed SA round trip differs");
            if (width > 32 && width < 64 &&
                packed.packedBytes().size() >= values.size()*sizeof(std::uint64_t))
                throw std::runtime_error("wide SA was stored as 64-bit entries");
            auto writable=PackedSuffixArray::zeroed(values.size(),width);
            auto lazy=PackedSuffixArray::lazyZeroed(values.size(),width);
            for (std::size_t i=values.size();i-->0;) lazy.set(i,values[i]);
            auto lazy_copy=lazy;
            lazy.set(0,0);
            auto lazy_move=std::move(lazy);
            if (!std::ranges::equal(lazy_copy,values) || lazy_move[0]!=0 ||
                lazy_move[values.size()-1]!=values.back())
                throw std::runtime_error("lazy packed ownership or copy independence differs");
            for (std::size_t i=values.size();i-->0;) writable.set(i,values[i]);
            // Overwriting must clear old bits without modifying adjacent ranks.
            for (std::size_t i=0;i<values.size();i+=2) writable.set(i,0);
            for (std::size_t i=0;i<values.size();++i)
                if (writable[i]!=(i%2 ? values[i] : 0))
                    throw std::runtime_error("random packed writes corrupted neighboring ranks");
            const auto encoded = packed.packedBytes();
            auto owner=std::make_shared<std::vector<std::uint8_t>>(encoded.begin(),encoded.end());
            owner->resize(owner->size()+7,0);
            auto view=PackedSuffixArray::fromView(width,values.size(),*owner,owner);
            auto copy=view;
            owner.reset(); view=PackedSuffixArray{};
            if (!std::ranges::equal(copy,values)) throw std::runtime_error("mapped owner lifetime lost");
            const auto reloaded = PackedSuffixArray::fromBytes(width, values.size(),
                std::vector<std::uint8_t>(encoded.begin(), encoded.end()));
            if (!std::ranges::equal(reloaded, values) ||
                reloaded[0] != values[0] ||
                reloaded[values.size()-1] != values.back())
                throw std::runtime_error("packed SA reload differs");
            // 64-aligned ranges written in reverse order, as parallel chunks may
            // finish, must produce the same bytes as one sequential pass.
            auto ranged=PackedSuffixArray::lazyZeroed(values.size(),width);
            for (std::size_t begin=(values.size()-1)/64*64;;begin-=64) {
                auto writer=ranged.rangeWriter(begin);
                for (std::size_t i=begin;i<std::min(values.size(),begin+64);++i)
                    writer.put(values[i]);
                writer.finish();
                if (begin==0) break;
            }
            if (!std::ranges::equal(ranged.packedBytes(), expected))
                throw std::runtime_error("range writers differ from sequential packing");
            try {
                ranged.rangeWriter(1);
                throw std::runtime_error("unaligned range writer accepted");
            } catch (const std::invalid_argument&) {}
        }
        const std::uint64_t human_physical = 6'200'000'000ULL;
        // In-place conversion must retain its original base address, overwrite
        // nonzero input correctly, preserve guards and survive owner copies.
        for (unsigned width : {1,7,8,9,31,32,33,34,35,39,40}) {
            const auto domain=std::uint64_t{1}<<width;
            for (std::size_t n : {0U,1U,2U,7U,63U,64U,65U,4095U,4096U,4097U,20001U}) {
                nustar::star2711b::detail::CompactSuffixArray source(n+1000);
                std::vector<std::uint64_t> expected(n);
                for (std::size_t i=0;i<source.size();++i) source[i]=nustar::star2711b::detail::Position40::mask;
                for (std::size_t i=0;i<n;++i) source[i]=expected[i]=(i*4294967291ULL+12345)%domain;
                source.resize(n);
                auto* original=reinterpret_cast<const std::uint8_t*>(source.begin());
                auto packed=PackedSuffixArray::packAndRelease(std::move(source),domain);
                if (packed.packedBytes().data()!=original || source.size()!=0)
                    throw std::runtime_error("compact packing did not transfer its allocation");
                if (!std::ranges::equal(packed,expected))
                    throw std::runtime_error("in-place packing overwrote unread positions");
                PackedSuffixArray reference(expected,domain);
                if (!std::ranges::equal(packed.packedBytes(),reference.packedBytes()))
                    throw std::runtime_error("in-place packed bytes or padding differ");
                for (unsigned j=0;j<7;++j) if (original[packed.packedBytes().size()+j])
                    throw std::runtime_error("in-place reader guards are not zero");
                auto copy=packed;
                packed=PackedSuffixArray{};
                if (!std::ranges::equal(copy,expected))
                    throw std::runtime_error("in-place packed owner lifetime failed");
            }
        }
        {
            // Cross the source-page release boundary, with values above 32 bits.
            std::vector<std::uint64_t> source((1U<<20)+257);
            for (std::size_t i=0;i<source.size();++i) source[i]=(i*104729ULL)%human_physical;
            const auto count=source.size();
            auto consumed=PackedSuffixArray::packAndRelease(std::move(source),human_physical);
            for (std::size_t i=0;i<count;++i)
                if (consumed[i]!=(i*104729ULL)%human_physical)
                    throw std::runtime_error("consuming pack lost values at a released source page");
        }
        const std::vector<std::uint64_t> human_values{
            0, 4'294'967'296ULL, human_physical-1, 4'900'000'001ULL};
        const PackedSuffixArray human(human_values, human_physical);
        if (human.bitWidth() != 33 ||
            human.packedBytes().size() != 17 ||
            !std::ranges::equal(human, human_values))
            throw std::runtime_error("human-scale 33-bit positions differ");
        const std::vector<std::uint32_t> legacy{0, 42, 1000};
        if (!std::ranges::equal(PackedSuffixArray(legacy, 1024), legacy))
            throw std::runtime_error("32-bit builder positions differ");
        const PackedSuffixArray wide_domain_from_32_bit_builder(legacy, human_physical);
        if (wide_domain_from_32_bit_builder.bitWidth() != 33 ||
            wide_domain_from_32_bit_builder.packedBytes().size() != 13 ||
            !std::ranges::equal(wide_domain_from_32_bit_builder, legacy))
            throw std::runtime_error("wide domain fell back to 64-bit entries");
        try {
            PackedSuffixArray::fromBytes(31, 2, {});
            throw std::runtime_error("truncated packed SA accepted");
        } catch (const std::invalid_argument&) {}
        if (PackedSuffixArray::fromBytes(3, 1, {0xff})[0] != 7)
            throw std::runtime_error("unused packed padding affected a value");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
