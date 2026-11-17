#pragma once
// Private ordered BGZF output. Shares the sorter's workers; no nested pool.
#include "io_util.hpp"
#include "work_queue.hpp"
#include <atomic>
#include <chrono>
#include <fcntl.h>

namespace nustar::star2711b::bam {
struct BgzfMetrics {
    IoCounters output;
    double prepare_seconds = 0, compression_seconds = 0, write_seconds = 0, wait_seconds = 0, sync_seconds = 0;
};
class BgzfWriter {
    using Clock = std::chrono::steady_clock;
    static double elapsed(Clock::time_point t) { return std::chrono::duration<double>(Clock::now()-t).count(); }
    struct Compressed { Buffer data; double prepare = 0, compress = 0; };
    static void compress(std::span<const std::byte> input, Buffer& output, int level);
    WorkQueue& pool_;
    std::atomic<bool>& stopped_;
    IoCalls calls_;
    int fd_ = -1, level_;
    std::size_t limit_, write_target_, ready_bytes_ = 0;
    std::deque<std::future<Compressed>> pending_;
    std::vector<Buffer> ready_;
    BgzfMetrics metrics_;
    void flushReady() {
        if (!ready_bytes_) return;
        std::vector<iovec> vectors;
        for (auto& b : ready_) vectors.push_back({b.data.get(), b.size});
        const auto start = Clock::now();
        vectorIo(fd_, std::move(vectors), 0, true, false, metrics_.output, calls_);
        metrics_.write_seconds += elapsed(start);
        ready_.clear(); ready_bytes_ = 0;
    }
  public:
    BgzfWriter(const std::string& path, int level, std::uint64_t budget, WorkQueue& pool,
               std::atomic<bool>& stopped, IoCalls calls = {})
        : pool_(pool), stopped_(stopped), calls_(calls), level_(level),
          limit_(std::max<std::uint64_t>(1, std::min<std::uint64_t>(2*pool.concurrency(), budget/(8ULL<<20)))),
          write_target_(std::max<std::uint64_t>(65536, std::min<std::uint64_t>(1ULL<<20,budget/32))) {
        fd_ = ::open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC, 0666);
        if (fd_ < 0) throw std::runtime_error("cannot open sorted BAM output");
    }
    ~BgzfWriter() {
        // Tasks may retain input arenas. Join before destroying their owner.
        for (auto& f : pending_) f.wait();
        if (fd_ >= 0) ::close(fd_);
    }
    template<class Owner, class RecordAt>
    void appendChunk(std::shared_ptr<Owner> owner, std::size_t record, std::size_t offset,
                     std::size_t length, RecordAt get) {
        if (!length) return;
        if (stopped_.load(std::memory_order_relaxed)) throw std::runtime_error("BAM output cancelled");
        while (pending_.size() >= limit_) drainOne();
        const auto level = level_;
        auto* stopped = &stopped_;
        auto future = pool_.submit([owner=std::move(owner), record, offset, length, get, level, stopped]() mutable {
            Compressed result;
            std::array<std::byte, 0xff00> input;
            std::size_t used = 0;
            auto start = Clock::now();
            while (length) {
                if (stopped->load(std::memory_order_relaxed)) throw std::runtime_error("BAM compression cancelled");
                const auto bytes = get(*owner, record);
                const auto n = std::min({length, bytes.size()-offset, input.size()-used});
                if (!n) throw std::runtime_error("invalid BAM output chunk");
                std::memcpy(input.data()+used, bytes.data()+offset,n);
                used += n; offset += n; length -= n;
                if (offset == bytes.size()) { ++record; offset = 0; }
                if (used == input.size() || !length) {
                    result.prepare += elapsed(start); start = Clock::now();
                    compress(std::span(input).first(used), result.data, level);
                    result.compress += elapsed(start); start = Clock::now(); used = 0;
                }
            }
            owner.reset(); // Release arena credits before publishing completion.
            return result;
        }, true);
        try { pending_.push_back(std::move(future)); }
        catch (...) { if (future.valid()) future.wait(); throw; }
    }
    void append(std::span<const std::byte> bytes) {
        while (!bytes.empty()) {
            const auto n = std::min<std::size_t>(bytes.size(),1ULL<<20);
            auto owned = std::make_shared<Buffer>(); owned->append(bytes.first(n));
            appendChunk(owned,0,0,n,[](const Buffer& b, std::size_t) { return b.bytes(); });
            bytes = bytes.subspan(n);
        }
    }
    bool drainOne() {
        if (pending_.empty()) return false;
        const auto start = Clock::now();
        auto future = std::move(pending_.front()); pending_.pop_front();
        auto result = future.get();
        metrics_.wait_seconds += elapsed(start);
        metrics_.prepare_seconds += result.prepare; metrics_.compression_seconds += result.compress;
        ready_bytes_ += result.data.size; ready_.push_back(std::move(result.data));
        if (ready_bytes_ >= write_target_) flushReady();
        return true;
    }
    void finish() {
        while (drainOne()) {}
        flushReady();
        // Exactly one standard BGZF EOF marker, independent of compression level.
        const unsigned char eof[28] = {31,139,8,4,0,0,0,0,0,255,6,0,66,67,2,0,27,0,3,0,0,0,0,0,0,0,0,0};
        const auto write_start = Clock::now();
        vectorIo(fd_, {{const_cast<unsigned char*>(eof),sizeof(eof)}},0,true,false,metrics_.output,calls_);
        metrics_.write_seconds += elapsed(write_start);
        const auto start = Clock::now();
        int code;
        do { code = calls_.sync(fd_); } while (code < 0 && errno == EINTR);
        metrics_.sync_seconds += elapsed(start);
        if (code < 0) throw std::runtime_error("cannot synchronize sorted BAM output");
        const auto fd = fd_; fd_ = -1;
        if (::close(fd) != 0) throw std::runtime_error("cannot close sorted BAM output");
    }
    const BgzfMetrics& metrics() const { return metrics_; }
};
}
