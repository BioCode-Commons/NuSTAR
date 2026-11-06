#pragma once

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <bit>
#include <cstring>
#include <iterator>
#include <memory>
#include <span>
#include <vector>

namespace nustar::star2711b {

// LSB-first, little-endian storage for SA positions and multilevel SAi values.
// In the suffix array, positions remain physical-sequence
// coordinates; unlike STAR's SA, their strand is already encoded by the
// physical coordinate and does not need a separate bit.
class PackedSuffixArray {
public:
    class const_iterator {
    public:
        using iterator_concept = std::forward_iterator_tag;
        using iterator_category = std::forward_iterator_tag;
        using value_type = std::uint64_t;
        using difference_type = std::ptrdiff_t;

        const_iterator() = default;
        value_type operator*() const noexcept { return (*owner_)[index_]; }
        const_iterator& operator++() noexcept { ++index_; return *this; }
        const_iterator operator++(int) noexcept {
            auto previous = *this;
            ++*this;
            return previous;
        }
        friend bool operator==(const const_iterator&, const const_iterator&) = default;

    private:
        friend class PackedSuffixArray;
        const_iterator(const PackedSuffixArray* owner, std::size_t index)
            : owner_(owner), index_(index) {}
        const PackedSuffixArray* owner_ = nullptr;
        std::size_t index_ = 0;
    };

    PackedSuffixArray() = default;
    PackedSuffixArray(const PackedSuffixArray&);
    PackedSuffixArray& operator=(const PackedSuffixArray&);
    PackedSuffixArray(PackedSuffixArray&&) noexcept = default;
    PackedSuffixArray& operator=(PackedSuffixArray&&) noexcept = default;
    PackedSuffixArray(std::span<const std::uint32_t> values,
                      std::uint64_t physical_length);
    PackedSuffixArray(std::span<const std::uint64_t> values,
                      std::uint64_t physical_length);
    static PackedSuffixArray fromBytes(std::uint8_t width, std::size_t count,
                                       std::vector<std::uint8_t> bytes);
    // Read-only owner-backed view; the caller supplies seven readable guard
    // bytes after the packed payload. Shared ownership pins the mapped file.
    static PackedSuffixArray fromView(std::uint8_t width, std::size_t count,
        std::span<const std::uint8_t> bytes, std::shared_ptr<const void> owner,
        bool file_mapping=false);
    static std::uint8_t widthFor(std::uint64_t physical_length);
    static std::size_t requiredBytes(std::size_t count, std::uint8_t width);

    // Construction-only write to owned storage with an in-range rank/value.
    // Mapped views are immutable; the multilevel prefix builder uses owned storage.
    static PackedSuffixArray zeroed(std::size_t count, std::uint8_t width);
    void set(std::size_t rank, std::uint64_t value);
    // calloc keeps untouched large pages uncommitted. Used by the streaming
    // junction merge, where an eagerly cleared second SA doubles peak RSS.
    static PackedSuffixArray lazyZeroed(std::size_t count,std::uint8_t width);
    // Consume construction-only positions. On Linux, return fully consumed
    // source pages while filling the packed array, instead of retaining both.
    static PackedSuffixArray packAndRelease(std::vector<std::uint64_t> values,
                                             std::uint64_t physical_length);
    // Advisory only: immutable file contents remain readable after discard.
    void discardReadPages(std::size_t first_rank,std::size_t past_rank) const;

    // Writes consecutive ranks of zeroed owned storage, one whole 64-bit word
    // at a time. The first rank must be a multiple of 64, so the range starts
    // on a word boundary; a range that also ends at a multiple of 64 (or at
    // the last rank) touches no byte outside itself. Writers of disjoint such
    // ranges can therefore run on different threads.
    class RangeWriter {
    public:
        void put(std::uint64_t value) noexcept {
            word_ |= value << bits_;
            bits_ += width_;
            if (bits_ >= 64) {
                store(word_);
                bits_ -= 64;
                word_ = bits_ ? value >> (width_-bits_) : 0;
            }
        }
        // Stores the bytes of a final partial word.
        void finish() noexcept {
            for (unsigned i = 0; 8*i < bits_; ++i)
                out_[i] = static_cast<std::uint8_t>(word_ >> (8*i));
            bits_ = 0; word_ = 0;
        }
    private:
        friend class PackedSuffixArray;
        RangeWriter(std::uint8_t* out, std::uint8_t width) : out_(out), width_(width) {}
        void store(std::uint64_t word) noexcept {
            if constexpr (std::endian::native == std::endian::little)
                std::memcpy(out_, &word, sizeof(word));
            else
                for (unsigned i = 0; i < 8; ++i)
                    out_[i] = static_cast<std::uint8_t>(word >> (8*i));
            out_ += 8;
        }
        std::uint8_t* out_;
        unsigned width_;
        unsigned bits_ = 0;
        std::uint64_t word_ = 0;
    };
    RangeWriter rangeWriter(std::size_t first_rank);

    // As in STAR PackedArray.h, expose this tiny hot accessor to the compiler
    // and retain the mask with the layout. memcpy permits unaligned loads
    // without violating C++ aliasing; widths above 57 keep the second-word path.
    std::uint64_t operator[](std::size_t rank) const noexcept {
        const auto read_word = [](const std::uint8_t* bytes) {
            std::uint64_t word = 0;
            if constexpr (std::endian::native == std::endian::little)
                std::memcpy(&word, bytes, sizeof(word));
            else
                for (unsigned i = 0; i < 8; ++i)
                    word |= std::uint64_t(bytes[i]) << (8*i);
            return word;
        };
        const auto bit = rank*width_;
        const auto shift = static_cast<unsigned>(bit%8);
        const auto* source = data_+bit/8;
        auto value = read_word(source) >> shift;
        if (static_cast<unsigned>(width_)+shift > 64)
            value |= read_word(source+8) << (64-shift);
        return value & mask_;
    }
    std::size_t size() const noexcept { return count_; }
    bool empty() const noexcept { return count_ == 0; }
    std::uint8_t bitWidth() const noexcept { return width_; }
    std::size_t storageBytes() const noexcept { return data_ ? payload_bytes_+7 : 0; }
    std::span<const std::uint8_t> packedBytes() const noexcept {
        return {data_, payload_bytes_};
    }
    const_iterator begin() const noexcept { return {this, 0}; }
    const_iterator end() const noexcept { return {this, count_}; }

private:
    PackedSuffixArray(std::size_t count, std::uint64_t physical_length);
    void put(std::size_t rank, std::uint64_t value);
    std::size_t count_ = 0;
    std::size_t payload_bytes_ = 0;
    std::uint8_t width_ = 1;
    std::uint64_t mask_ = 1;
    // Seven zero bytes make unaligned 64-bit reads safe for the final entry.
    std::vector<std::uint8_t> bytes_;
    std::unique_ptr<std::uint8_t,decltype(&std::free)> lazy_bytes_{nullptr,&std::free};
    bool file_mapping_=false;
    std::shared_ptr<const void> mapped_owner_;
    const std::uint8_t* data_ = nullptr;
};

} // namespace nustar::star2711b
