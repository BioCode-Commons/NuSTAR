#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>
#include "input_config.hpp"

namespace nustar::star2711b {

enum class FastxFormat : std::uint8_t { Fasta, Fastq };
struct RawMate {
    std::string header, sequence;
    std::optional<std::string> quality; // Null for FASTA; original byte values for FASTQ.
    std::string header_extra;
};
struct RawFragment {
    FragmentId id;
    Layout layout;
    FastxFormat format; // Detected per fragment from mate 1, matching STAR's readFileType.
    std::array<std::optional<RawMate>, 2> mates; // Presence is explicit; empty sequence is valid.
    char legacy_read_filter;
    std::optional<RawMate> barcode;
};
struct ReadBatch {
    BatchId id;
    std::vector<RawFragment> fragments; // Owns complete paired records in input order.
};
class IReadSource {
public:
    virtual ~IReadSource() = default;
    // Called before workers start. A capable source prepares concurrent next()
    // calls; input bytes are read on demand by those callers, without an input mutex.
    virtual Status prepareParallel(std::uint64_t, std::uint32_t,
                                   std::uint64_t = UINT64_MAX) { return std::monostate{}; }
    virtual bool concurrentReads() const noexcept { return false; }
    // Null batch means clean EOF; an unmatched mate EOF is a Failure.
    virtual Result<std::optional<ReadBatch>> next(std::uint64_t target_bytes) = 0;
    virtual void cancel() noexcept = 0;
};
// Backends own only byte transport and FASTQ framing. A different IO strategy
// can be installed without changing the read preparation or mapper.
class IReadSourceBackend {
public:
    virtual ~IReadSourceBackend() = default;
    virtual Result<std::unique_ptr<IReadSource>> open(const InputSpec&) const = 0;
};
const IReadSourceBackend& fastqSourceBackend();
// Compatibility aliases for the same FASTQ backend.
const IReadSourceBackend& memoryMappedFastqBackend();
const IReadSourceBackend& bufferedFastqBackend();
Result<std::unique_ptr<IReadSource>> openIndexedFastqSource(
    const InputSpec&, std::uint64_t target_bytes, std::uint64_t limit);
Result<std::unique_ptr<IReadSource>> openReadSource(const InputSpec&);

// RawFragment is a transport envelope only; typed batches enter alignment preparation.
struct RawSingleRead {
    FragmentId id;
    FastxFormat format;
    RawMate mate;
    char legacy_read_filter;
};
struct RawReadPair {
    FragmentId id;
    FastxFormat format;
    std::array<RawMate, 2> mates;
    char legacy_read_filter;
};
struct SingleReadBatch { BatchId id; std::vector<RawSingleRead> reads; };
struct PairedReadBatch { BatchId id; std::vector<RawReadPair> pairs; };
Result<SingleReadBatch> takeSingleBatch(ReadBatch&&); // Reject a PE envelope before preparation.
Result<PairedReadBatch> takePairedBatch(ReadBatch&&); // Reject incomplete or SE envelopes.

} // namespace nustar::star2711b
