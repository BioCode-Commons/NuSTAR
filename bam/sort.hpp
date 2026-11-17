#pragma once
// Private implementation of coordinate sorting. Not part of bam/io.hpp.
#include "codec.hpp"
#include "io_util.hpp"
#include "../output/sink.hpp"
#include <functional>
#include <memory>

namespace nustar::star2711b::bam {
struct SortMetrics {
    std::uint64_t records = 0, spilled_bytes = 0, read_bytes = 0, partitions = 0;
    double radix_seconds = 0, pack_seconds = 0, output_seconds = 0;
    IoCounters temporary_reads, temporary_writes, output_writes;
    double read_seconds = 0, parse_seconds = 0, write_seconds = 0;
    double compression_seconds = 0, wait_seconds = 0, sync_seconds = 0;
};
class CoordinateSorter {
    struct Impl;
    std::unique_ptr<Impl> impl_;
  public:
    CoordinateSorter(std::filesystem::path directory, std::uint64_t budget, unsigned threads);
    ~CoordinateSorter();
    // Concurrent calls own independent slabs; no caller buffer is retained.
    void submit(const OutputBatch& batch);
    void finish();
    void emit(const std::function<void(std::uint64_t, std::span<const std::byte>)>& output,
              const std::function<void(std::span<const std::byte>)>& blocks = {});
    // Module-private output path; the public sink/config contracts are unchanged.
    void writeBam(const std::string& path, std::span<const std::byte> header, int compression,
                  const std::function<std::pair<std::string,std::string>(std::uint64_t)>& correction = {});
    void cancel() noexcept;
    SortMetrics metrics() const;
};
} // namespace nustar::star2711b::bam
