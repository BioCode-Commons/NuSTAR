#pragma once
#include <algorithm>
#include <array>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <system_error>
#include <sys/uio.h>
#include <unistd.h>
#include <vector>

namespace nustar::star2711b::bam {
// Reusable storage: growing does not clear bytes which the caller overwrites.
struct Buffer {
    std::unique_ptr<std::byte[]> data;
    std::size_t size = 0, capacity = 0;
    void reserve(std::size_t n) {
        if (n <= capacity) return;
        auto next = std::make_unique_for_overwrite<std::byte[]>(n);
        if (size) std::memcpy(next.get(), data.get(), size);
        data = std::move(next); capacity = n;
    }
    void resize(std::size_t n) { reserve(n); size = n; }
    void append(std::span<const std::byte> bytes) {
        if (bytes.size() > std::numeric_limits<std::size_t>::max() - size)
            throw std::runtime_error("BAM buffer size overflow");
        const auto n = size + bytes.size();
        if (n > capacity) reserve(std::max(n, capacity <= SIZE_MAX / 2 ? 2 * capacity : n));
        if (!bytes.empty()) std::memcpy(data.get() + size, bytes.data(), bytes.size());
        size = n;
    }
    std::span<std::byte> bytes() { return {data.get(), size}; }
    std::span<const std::byte> bytes() const { return {data.get(), size}; }
};
struct IoCounters {
    std::uint64_t calls = 0, bytes = 0;
    // <=4K,16K,64K,256K,1M,4M,16M,larger (actual syscall requests).
    std::array<std::uint64_t, 8> requests{};
    void request(std::size_t n) {
        unsigned bin = 0; std::size_t end = 4096;
        while (bin < 7 && n > end) { ++bin; end *= 4; }
        ++requests[bin]; ++calls;
    }
    void add(const IoCounters& other) {
        calls += other.calls; bytes += other.bytes;
        for (unsigned i = 0; i < requests.size(); ++i) requests[i] += other.requests[i];
    }
};
// Per-instance syscall injection is used only by module tests.
struct IoCalls {
    decltype(&::preadv) read = ::preadv;
    decltype(&::pwritev) pwrite = ::pwritev;
    decltype(&::writev) write = ::writev;
    decltype(&::fsync) sync = ::fsync;
};
inline void vectorIo(int fd, std::vector<iovec> vectors, std::uint64_t offset,
                     bool writing, bool positioned, IoCounters& stats, const IoCalls& calls = {}) {
    const auto maximum = std::uint64_t(std::numeric_limits<off_t>::max());
    std::uint64_t total = 0;
    for (const auto& v : vectors) {
        if (v.iov_len > maximum - total) throw std::runtime_error("BAM I/O length overflow");
        total += v.iov_len;
    }
    if (positioned && (offset > maximum || total > maximum - offset))
        throw std::runtime_error("BAM I/O offset overflow");
    const long configured = ::sysconf(_SC_IOV_MAX);
    const std::size_t limit = configured > 0 ? std::min<long>(configured, 1024) : 16;
    std::size_t first = 0;
    while (first < vectors.size()) {
        if (!vectors[first].iov_len) { ++first; continue; }
        std::vector<iovec> request;
        std::size_t bytes = 0;
        for (auto i = first; i < vectors.size() && request.size() < limit && bytes < (1ULL << 30); ++i) {
            auto v = vectors[i]; v.iov_len = std::min<std::size_t>(v.iov_len, (1ULL << 30) - bytes);
            request.push_back(v); bytes += v.iov_len;
        }
        stats.request(bytes);
        const auto result = writing ? (positioned ? calls.pwrite(fd, request.data(), request.size(), offset)
                                                 : calls.write(fd, request.data(), request.size()))
                                    : calls.read(fd, request.data(), request.size(), offset);
        if (result < 0 && errno == EINTR) continue;
        if (result < 0) {
            const int error = errno;
            throw std::system_error(error, std::generic_category(),
                std::string(writing ? (positioned ? "BAM temporary pwritev" : "BAM output writev") : "BAM temporary preadv") +
                " at offset " + std::to_string(offset) + " requesting " + std::to_string(bytes) + " bytes");
        }
        if (!result) throw std::runtime_error(writing ? "zero-length BAM write" : "truncated BAM sort data");
        if (std::uint64_t(result) > bytes) throw std::runtime_error("invalid BAM I/O result");
        stats.bytes += result; offset += result;
        std::size_t remaining = result;
        while (first < vectors.size() && remaining >= vectors[first].iov_len) {
            remaining -= vectors[first].iov_len; ++first;
        }
        if (remaining) {
            auto& v = vectors[first]; v.iov_base = static_cast<std::byte*>(v.iov_base) + remaining;
            v.iov_len -= remaining;
        }
    }
}
}
