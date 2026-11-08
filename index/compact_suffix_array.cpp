#include "compact_suffix_array.hpp"
#include <limits>
#include <new>
#include <stdexcept>
#include <cstdlib>
#if defined(__unix__)
#include <sys/mman.h>
#include <unistd.h>
#endif
namespace nustar::star2711b::detail {
CompactSuffixArray::CompactSuffixArray(std::size_t n):count_(n) {
    if (n>Position40::flag-1 || n>(std::numeric_limits<std::size_t>::max()-7)/5)
        throw std::length_error("compact suffix workspace exceeds its capacity");
    bytes_=n*5+7;
    if (!bytes_) return;
#if defined(__unix__)
    auto* p=::mmap(nullptr,bytes_,PROT_READ|PROT_WRITE,MAP_PRIVATE|MAP_ANONYMOUS,-1,0);
    if(p==MAP_FAILED) throw std::bad_alloc();
    data_=static_cast<Position40*>(p);
#else
    data_=static_cast<Position40*>(std::calloc(bytes_,1));
    if(!data_) throw std::bad_alloc();
#endif
}
void CompactSuffixArray::release() noexcept {
    if(!data_) return;
#if defined(__unix__)
    ::munmap(data_,bytes_);
#else
    std::free(data_);
#endif
    data_=nullptr;count_=bytes_=released_=0;
}
CompactSuffixArray::~CompactSuffixArray(){release();}
CompactSuffixArray::CompactSuffixArray(CompactSuffixArray&& o) noexcept {
    data_=std::exchange(o.data_,nullptr);count_=std::exchange(o.count_,0);
    bytes_=std::exchange(o.bytes_,0);released_=std::exchange(o.released_,0);
}
CompactSuffixArray& CompactSuffixArray::operator=(CompactSuffixArray&& o) noexcept {
    if(this!=&o){release();data_=std::exchange(o.data_,nullptr);count_=std::exchange(o.count_,0);
    bytes_=std::exchange(o.bytes_,0);released_=std::exchange(o.released_,0);}return *this;
}
void CompactSuffixArray::discardBefore(std::size_t n) noexcept {
#if defined(__unix__)
    const auto page=std::size_t(::sysconf(_SC_PAGESIZE));
    const auto end=(n*5/page)*page;
    if(end>released_){::madvise(reinterpret_cast<char*>(data_)+released_,end-released_,MADV_DONTNEED);released_=end;}
#else
    (void)n;
#endif
}

std::shared_ptr<const void> CompactSuffixArray::takePackedStorage(std::size_t retained_bytes) {
    if (!data_ || retained_bytes > bytes_ || released_ != 0)
        throw std::invalid_argument("packed storage must fit an undiscarded workspace");
    struct Owner {
        void* address=nullptr;
        std::size_t bytes=0;
        ~Owner() {
            if (!address) return;
#if defined(__unix__)
            ::munmap(address,bytes);
#else
            std::free(address);
#endif
        }
    };
    // Allocate the control block before transferring ownership (exception safe).
    auto owner=std::make_shared<Owner>();
#if defined(__unix__)
    const auto page=std::size_t(::sysconf(_SC_PAGESIZE));
    const auto retained_pages=(retained_bytes+page-1)/page*page;
    if (retained_pages < bytes_ &&
        ::munmap(reinterpret_cast<char*>(data_)+retained_pages,bytes_-retained_pages)==0)
        bytes_=retained_pages;
#endif
    owner->address=std::exchange(data_,nullptr);
    owner->bytes=std::exchange(bytes_,0);
    count_=released_=0;
    return std::shared_ptr<const void>(owner,owner->address);
}
}
