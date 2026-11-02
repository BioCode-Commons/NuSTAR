#pragma once

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace nustar::star2711b {

// Distinct tags prevent accidental interchange of coordinates and identities.
template<class Tag> struct U64 { std::uint64_t value; };
using GenomeOffset = U64<struct GenomeOffsetTag>;
using ReadOffset = U64<struct ReadOffsetTag>;
using BinId = U64<struct BinIdTag>;
using ChromosomeId = U64<struct ChromosomeIdTag>;
using GeneId = U64<struct GeneIdTag>;
using TranscriptId = U64<struct TranscriptIdTag>;
using JunctionId = U64<struct JunctionIdTag>;
using CandidateId = U64<struct CandidateIdTag>;
using WindowId = U64<struct WindowIdTag>;
using BatchId = U64<struct BatchIdTag>;
using ReferenceGeneration = U64<struct ReferenceGenerationTag>;
using IndexGeneration = U64<struct IndexGenerationTag>;
using OccurrenceSetId = U64<struct OccurrenceSetIdTag>;
using OccurrenceSourceId = U64<struct OccurrenceSourceIdTag>;
using Score = std::int32_t; // STAR intScore; do not widen floating-point expressions.

struct FragmentId {
    std::uint32_t file_group;
    std::uint64_t ordinal; // Zero-based position across all input groups; STAR iReadAll = ordinal + 1.
};

struct GenomeRange { GenomeOffset begin, end; };
struct ReadRange { ReadOffset begin, end; };
enum class MateId : std::uint8_t { First = 0, Second = 1 };
enum class Layout : std::uint8_t { SE = 1, PE = 2 };
enum class Orientation : std::uint8_t { Forward = 0, Reverse = 1 };
enum class Strand : std::uint8_t { Unknown = 0, Plus = 1, Minus = 2 };
enum class Motif : std::uint8_t {
    Noncanonical = 0, GT_AG = 1, CT_AC = 2, GC_AG = 3,
    CT_GC = 4, AT_AC = 5, GT_AT = 6
};
enum class BaseCode : std::uint8_t {
    A = 0, C = 1, G = 2, T = 3, N = 4, GenomePadding = 5, MateSpacer = 11
};
enum class Stage {
    Plan, Reference, Gtf, Sjdb, Input, Prepare, Seed, Window,
    Stitch, Overlap, Select, Filter, Project, Commit, SjReduce, Finish
};
enum class ErrorKind { Input, Config, Index, Capacity, Io, Cancelled, Internal };
struct Failure {
    ErrorKind kind;
    Stage stage;
    std::optional<FragmentId> fragment;
    int legacy_exit_code; // Preserve the STAR category when a legacy equivalent exists.
    std::string message;
};
template<class T> using Result = std::variant<T, Failure>;
using Status = Result<std::monostate>;

struct Diagnostic {
    enum class Level { Warning, Info } level;
    Stage stage;
    std::filesystem::path source;
    std::optional<std::uint64_t> line;
    std::string code;
    std::string message;
};
struct FileIdentity {
    std::filesystem::path path;
    std::uint64_t bytes;
    std::string sha256;
};
struct BaselineIdentity {
    std::string release; // Must be 2.7.11b for this compatibility profile.
    std::string source_id;
    std::string binary_sha256;
    std::string compiler_and_flags;
    std::string standard_library;
    std::vector<FileIdentity> dependencies;
};
struct OptionProvenance {
    std::string name;
    std::vector<std::string> raw_tokens;
    std::vector<std::string> resolved_tokens;
    std::string origin; // Default, parameter file, CLI, index, or derived value.
};
enum class OptionHandling { Projected, LegacyPath, Rejected };
struct OptionDecision {
    std::string name; // One decision for every effective STAR option, including defaults.
    OptionHandling handling;
    std::string owner_or_reason; // Module for Projected; explicit reason otherwise.
};

} // namespace nustar::star2711b
