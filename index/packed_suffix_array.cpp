#include "packed_suffix_array.hpp"

#include <bit>
#include <cstring>
#include <limits>
#include <stdexcept>
#if defined(__unix__)
#include <sys/mman.h>
#include <unistd.h>
#endif

namespace nustar::star2711b {
namespace {

std::uint64_t readWord(const std::uint8_t* bytes) noexcept {
    if constexpr (std::endian::native == std::endian::little) {
        std::uint64_t word;
        std::memcpy(&word, bytes, sizeof(word));
        return word;
    } else {
        std::uint64_t word = 0;
        for (unsigned i = 0; i < 8; ++i) word |= std::uint64_t(bytes[i]) << (8*i);
        return word;
    }
}

void writeWord(std::uint8_t* bytes, std::uint64_t word) noexcept {
    if constexpr (std::endian::native == std::endian::little)
        std::memcpy(bytes, &word, sizeof(word));
    else
        for (unsigned i = 0; i < 8; ++i)
            bytes[i] = static_cast<std::uint8_t>(word >> (8*i));
}


} // namespace

PackedSuffixArray::PackedSuffixArray(const PackedSuffixArray& other)
    : count_(other.count_), payload_bytes_(other.payload_bytes_), width_(other.width_),
      mask_(other.mask_), bytes_(other.lazy_bytes_ ? std::vector<std::uint8_t>(other.data_,other.data_+other.payload_bytes_+7) : other.bytes_),
      file_mapping_(other.file_mapping_), mapped_owner_(other.mapped_owner_),
      data_(mapped_owner_ ? other.data_ : bytes_.data()) {}
PackedSuffixArray& PackedSuffixArray::operator=(const PackedSuffixArray& other) {
    if (this!=&other) *this=PackedSuffixArray(other);
    return *this;
}

PackedSuffixArray PackedSuffixArray::fromView(std::uint8_t width, std::size_t count,
    std::span<const std::uint8_t> bytes, std::shared_ptr<const void> owner,bool file_mapping) {
    const auto payload=requiredBytes(count,width);
    if (!owner || bytes.size()<payload+7)
        throw std::invalid_argument("mapped packed array lacks owner or guards");
    PackedSuffixArray result;
    result.count_=count; result.payload_bytes_=payload; result.width_=width;
    result.mask_=~std::uint64_t{0}>>(64-width);
    result.data_=bytes.data(); result.mapped_owner_=std::move(owner);
    result.file_mapping_=file_mapping;
    return result;
}

std::uint8_t PackedSuffixArray::widthFor(std::uint64_t physical_length) {
    if (physical_length == 0)
        throw std::invalid_argument("physical sequence length must be positive");
    const auto bits = std::bit_width(physical_length-1);
    return static_cast<std::uint8_t>(bits == 0 ? 1 : bits);
}

std::size_t PackedSuffixArray::requiredBytes(std::size_t count,
                                              std::uint8_t width) {
    if (width == 0 || width > 64 ||
        count > (std::numeric_limits<std::size_t>::max()-7)/width)
        throw std::invalid_argument("invalid packed suffix array dimensions");
    return (count*width+7)/8;
}

PackedSuffixArray::PackedSuffixArray(std::size_t count,
                                     std::uint64_t physical_length)
    : count_(count), payload_bytes_(requiredBytes(count, widthFor(physical_length))),
      width_(widthFor(physical_length)),
      mask_(std::numeric_limits<std::uint64_t>::max() >> (64-width_)),
      bytes_(payload_bytes_+7, 0) {
    data_=bytes_.data();
}

PackedSuffixArray::PackedSuffixArray(std::span<const std::uint32_t> values,
                                     std::uint64_t physical_length)
    : PackedSuffixArray(values.size(), physical_length) {
    for (std::size_t i = 0; i < values.size(); ++i)
        put(i, values[i]);
}

PackedSuffixArray::PackedSuffixArray(std::span<const std::uint64_t> values,
                                     std::uint64_t physical_length)
    : PackedSuffixArray(values.size(), physical_length) {
    for (std::size_t i = 0; i < values.size(); ++i)
        put(i, values[i]);
}

void PackedSuffixArray::put(std::size_t rank, std::uint64_t value) {
    const auto bit = rank*width_;
    const auto shift = static_cast<unsigned>(bit%8);
    auto* target = const_cast<std::uint8_t*>(data_)+bit/8;
    writeWord(target, readWord(target) | (value << shift));
    if (static_cast<unsigned>(width_)+shift > 64)
        writeWord(target+8, readWord(target+8) | (value >> (64-shift)));
}

PackedSuffixArray PackedSuffixArray::zeroed(std::size_t count, std::uint8_t width) {
    auto bytes = std::vector<std::uint8_t>(requiredBytes(count, width)+7, 0);
    bytes.resize(requiredBytes(count, width));
    return fromBytes(width, count, std::move(bytes));
}

void PackedSuffixArray::set(std::size_t rank, std::uint64_t value) {
    const auto bit = rank*width_;
    const auto shift = unsigned(bit%8);
    auto* target = const_cast<std::uint8_t*>(data_)+bit/8;
    writeWord(target, (readWord(target) & ~(mask_ << shift)) | (value << shift));
    if (unsigned(width_)+shift > 64)
        writeWord(target+8, (readWord(target+8) & ~(mask_ >> (64-shift))) |
            (value >> (64-shift)));
}

PackedSuffixArray PackedSuffixArray::lazyZeroed(std::size_t count,std::uint8_t width) {
    PackedSuffixArray result;
    result.count_=count;result.width_=width;
    result.payload_bytes_=requiredBytes(count,width);
    result.mask_=~std::uint64_t{0}>>(64-width);
    result.lazy_bytes_.reset(static_cast<std::uint8_t*>(std::calloc(result.payload_bytes_+7,1)));
    if(!result.lazy_bytes_) throw std::bad_alloc();
    result.data_=result.lazy_bytes_.get();
    return result;
}
void PackedSuffixArray::discardReadPages(std::size_t first,std::size_t past) const {
#if defined(__unix__)
    if(!file_mapping_) return;
    const auto page=static_cast<std::uintptr_t>(::sysconf(_SC_PAGESIZE));
    const auto base=reinterpret_cast<std::uintptr_t>(data_);
    const auto begin=(base+first*width_/8+page-1)/page*page;
    const auto end=(base+past*width_/8)/page*page;
    if(end>begin) ::madvise(reinterpret_cast<void*>(begin),end-begin,MADV_DONTNEED);
#else
    (void)first;(void)past;
#endif
}

PackedSuffixArray PackedSuffixArray::fromBytes(std::uint8_t width,
    std::size_t count, std::vector<std::uint8_t> bytes) {
    const auto payload = requiredBytes(count, width);
    if (bytes.size() != payload)
        throw std::invalid_argument("packed suffix array byte count differs");
    PackedSuffixArray result;
    result.count_ = count;
    result.width_ = width;
    result.mask_ = std::numeric_limits<std::uint64_t>::max() >> (64-width);
    result.payload_bytes_ = payload;
    result.bytes_ = std::move(bytes);
    result.bytes_.resize(payload+7, 0);
    result.data_=result.bytes_.data();
    return result;
}

} // namespace nustar::star2711b
