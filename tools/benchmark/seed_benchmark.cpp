// Reproducible seed-stage benchmark, independent of FASTQ parsing and tracing.
//
//   seed-benchmark INDEX FIRST.fastq -- 5 --reads 100000
//   seed-benchmark INDEX FIRST.fastq SECOND.fastq 5 --reads 100000
//
// Build this same source against both revisions. Alternate process order when
// comparing binaries, and use the same index, reads, compiler and build flags.
// Reads are encoded once, limited to --reads fragments (default 100000). Each
// round times --sweeps whole sweeps using the production storage-reuse search
// API and reports both their total time and the mean seconds per sweep.
// A separate untimed sweep hashes every output field with SHA-256 and checks
// it against an initial warmup sweep. Parsing, hashing, index loading and JSON
// output are excluded from seconds. Peak RSS includes all of those operations
// and the preloaded reads; it is a whole-process high-water mark, in bytes.
// This is a single-threaded microbenchmark, not an end-to-end mapping benchmark.

#include "../../seed_search/star_search.hpp"
#include "../../common/sha256.hpp"

#include <array>
#include <charconv>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

using namespace nustar::star2711b;

namespace {

struct FastqRecord { std::string name, sequence; };

bool readFastq(std::ifstream& input, FastqRecord& record) {
    std::string plus, quality;
    if (!std::getline(input, record.name)) {
        if (!input.eof()) throw std::runtime_error("FASTQ read failed");
        return false;
    }
    if (!std::getline(input, record.sequence) || !std::getline(input, plus) ||
        !std::getline(input, quality) || record.name.empty() ||
        record.name[0] != '@' || plus.empty() || plus[0] != '+' ||
        quality.size() != record.sequence.size())
        throw std::runtime_error("invalid FASTQ record");
    return true;
}

BaseCode encode(char base) {
    switch (base) {
    case 'A': case 'a': return BaseCode::A;
    case 'C': case 'c': return BaseCode::C;
    case 'G': case 'g': return BaseCode::G;
    case 'T': case 't': return BaseCode::T;
    default: return BaseCode::N;
    }
}

BaseCode complement(BaseCode base) {
    const auto code = static_cast<std::uint8_t>(base);
    return code < 4 ? static_cast<BaseCode>(code^3) : base;
}

// Keep the layout identical to seed-trace: mate 2 is reverse-complemented,
// separated by MateSpacer, and the remaining views preserve STAR orientation.
MappingRead prepare(const FastqRecord& first, const FastqRecord* second,
                    std::uint64_t ordinal) {
    MappingRead read{};
    read.fragment = FragmentId{0, ordinal};
    auto& sequence = read.encoded[0];
    sequence.reserve(first.sequence.size() + (second ? second->sequence.size()+1 : 0));
    for (char base : first.sequence) sequence.push_back(encode(base));
    read.mates.push_back(PreparedMate{
        first.sequence.size(), first.sequence.size(), {},
        ReadRange{ReadOffset{0}, ReadOffset{first.sequence.size()}}, false});
    if (second) {
        sequence.push_back(BaseCode::MateSpacer);
        for (auto it = second->sequence.rbegin(); it != second->sequence.rend(); ++it)
            sequence.push_back(complement(encode(*it)));
        const auto begin = first.sequence.size()+1;
        read.mates.push_back(PreparedMate{
            second->sequence.size(), second->sequence.size(), {},
            ReadRange{ReadOffset{begin}, ReadOffset{sequence.size()}}, true});
    }
    read.encoded[1].reserve(sequence.size());
    for (auto base : sequence) read.encoded[1].push_back(complement(base));
    read.encoded[2].assign(read.encoded[1].rbegin(), read.encoded[1].rend());
    read.budget = ReadBudget{0, sequence.size(), sequence.size(),
                            first.sequence.size()+(second ? second->sequence.size() : 0)};
    return read;
}

template<class T> T checked(Result<T> result) {
    if (const auto* failure = std::get_if<Failure>(&result))
        throw std::runtime_error(failure->message);
    return std::get<T>(std::move(result));
}

std::uint64_t number(std::string_view value, std::string_view label) {
    std::uint64_t result = 0;
    const auto parsed = std::from_chars(value.data(), value.data()+value.size(), result);
    if (parsed.ec != std::errc{} || parsed.ptr != value.data()+value.size())
        throw std::runtime_error("invalid unsigned integer for " + std::string(label));
    return result;
}

// Canonical little-endian uint64 fields, not raw object bytes (which may have
// padding). Vector sizes and per-read ordinals delimit every result exactly.
class Digest {
public:
    Digest() { hash_.update("NuSTAR-seed-result-v1"); }
    void integer(std::uint64_t value) {
        if (used_ == bytes_.size()) flush();
        for (unsigned byte = 0; byte < 8; ++byte)
            bytes_[used_++] = std::byte((value >> (8*byte)) & 255);
    }
    void result(std::uint64_t ordinal, const StarSeedSearchResult& result) {
        integer(ordinal);
        integer(result.map_marker);
        const auto& c = result.counters;
        integer(c.all_alignments); integer(c.unique_alignments);
        integer(c.multiple_alignments); integer(c.min_multimaps);
        integer(c.min_multimap_length); integer(c.max_unique_length);
        integer(c.no_good_piece_length);
        integer(result.splits.size());
        for (const auto& split : result.splits) {
            integer(split.read.begin.value); integer(split.read.end.value);
            integer(static_cast<unsigned>(split.mate));
        }
        integer(result.pieces.size());
        for (const auto& piece : result.pieces) {
            integer(piece.read_start); integer(piece.length);
            integer(piece.direction); integer(piece.multiplicity);
            integer(piece.sa_begin); integer(piece.sa_end); integer(piece.fragment);
        }
    }
    std::string finish() { flush(); return hash_.finish(); }
private:
    void flush() {
        hash_.update(std::span<const std::byte>(bytes_.data(), used_));
        used_ = 0;
    }
    Sha256 hash_;
    std::array<std::byte,4096> bytes_{};
    std::size_t used_ = 0;
};

struct Validation { std::string digest; std::uint64_t pieces; };

Validation validate(const StarSeedSearch& search, const std::vector<MappingRead>& reads,
                    StarSeedSearchResult& result) {
    Digest digest;
    digest.integer(reads.size());
    std::uint64_t pieces = 0;
    for (const auto& read : reads) {
        checked(search.search(read, result));
        digest.result(read.fragment.ordinal, result);
        pieces += result.pieces.size();
    }
    return {digest.finish(), pieces};
}

std::uint64_t peakRssBytes() {
#if defined(__unix__) || defined(__APPLE__)
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0)
        throw std::runtime_error("getrusage failed");
#if defined(__APPLE__)
    return static_cast<std::uint64_t>(usage.ru_maxrss);
#else
    return static_cast<std::uint64_t>(usage.ru_maxrss)*1024;
#endif
#else
    return 0;
#endif
}

void usage() {
    std::cerr << "Usage: seed-benchmark INDEX FIRST.fastq SECOND.fastq-or-- REPEATS\n"
        "  [--reads N] [--sweeps N] [--seedSearchLmax N] [--seedSearchStartLmax N]\n"
        "  [--seedMultimapNmax N] [--seedSplitMin N] [--seedMapMin N]\n"
        "Use -- in place of SECOND.fastq for single-end input. --reads defaults\n"
        "to 100000 fragments; --sweeps defaults to 1. REPEATS, --reads and\n"
        "--sweeps must be positive. Each round has bulk timed search sweeps\n"
        "and one untimed full-result digest sweep. seconds is time per sweep.\n";
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 5) { usage(); return 2; }
    try {
        const std::filesystem::path directory = argv[1];
        const bool paired = std::string_view(argv[3]) != "--";
        const auto repeats = number(argv[4], "REPEATS");
        std::uint64_t limit = 100000;
        std::uint64_t sweeps = 1;
        SeedConfig config{50, 0, 1.0, 10000, 1000, 12, 5, 10, 0};
        for (int i = 5; i < argc;) {
            const std::string option = argv[i++];
            if (i == argc) throw std::runtime_error("missing value for " + option);
            const auto value = number(argv[i++], option);
            if (option == "--reads") limit = value;
            else if (option == "--sweeps") sweeps = value;
            else if (option == "--seedSearchLmax") config.seedSearchLmax = value;
            else if (option == "--seedSearchStartLmax") config.seedSearchStartLmax = value;
            else if (option == "--seedMultimapNmax") config.seedMultimapNmax = value;
            else if (option == "--seedSplitMin") config.seedSplitMin = value;
            else if (option == "--seedMapMin") config.seedMapMin = value;
            else throw std::runtime_error("unsupported option: " + option);
        }
        if (!limit || !repeats || !sweeps || limit > std::numeric_limits<std::size_t>::max())
            throw std::runtime_error("REPEATS, --reads and --sweeps must be positive and fit this platform");

        std::vector<MappingRead> reads;
        std::uint64_t input_bases = 0;
        {
            std::ifstream first(argv[2]);
            std::ifstream second(paired ? argv[3] : "");
            if (!first || (paired && !second)) throw std::runtime_error("cannot open FASTQ");
            FastqRecord a, b;
            // Do not reserve the requested limit: a short input should only
            // allocate storage for the fragments it actually contains.
            while (reads.size() < limit && readFastq(first, a)) {
                if (paired && !readFastq(second, b))
                    throw std::runtime_error("FASTQ mates differ in record count");
                input_bases += a.sequence.size() + (paired ? b.sequence.size() : 0);
                reads.push_back(prepare(a, paired ? &b : nullptr, reads.size()));
            }
            if (reads.size() < limit && paired && readFastq(second, b))
                throw std::runtime_error("FASTQ mates differ in record count");
        }
        if (reads.empty()) throw std::runtime_error("FASTQ contains no reads");
        auto reference = checked(loadGenomeReference(directory));
        SuffixArrayIndexBackend backend;
        SearchIndexRequest request{"nustar-sa", directory/"index.sa", {}};
        auto index = std::dynamic_pointer_cast<const SuffixArrayIndex>(
            checked(backend.load(request, reference)));
        if (!index) throw std::runtime_error("index type mismatch");
        StarSeedSearch search(index, reference, config);
        StarSeedSearchResult result;
        const auto expected = validate(search, reads, result); // Untimed warmup.
        struct Round { double seconds, total_seconds, cpu_seconds; std::string digest; };
        std::vector<Round> rounds;
        for (std::uint64_t round = 0; round < repeats; ++round) {
            const auto cpu_started = std::clock();
            const auto started = std::chrono::steady_clock::now();
            for (std::uint64_t sweep = 0; sweep < sweeps; ++sweep) {
                std::uint64_t pieces = 0;
                for (const auto& read : reads) {
                    checked(search.search(read, result));
                    pieces += result.pieces.size();
                }
                // Check each sweep separately, avoiding expected.pieces*sweeps
                // overflow and retaining an observable result from every pass.
                if (pieces != expected.pieces)
                    throw std::runtime_error("seed piece count changed between timed sweeps");
            }
            const double total_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now()-started).count();
            const auto cpu_finished = std::clock();
            if (cpu_started==std::clock_t(-1) || cpu_finished==std::clock_t(-1))
                throw std::runtime_error("process CPU clock unavailable");
            const double cpu_seconds = double(cpu_finished-cpu_started)/CLOCKS_PER_SEC/sweeps;
            const auto actual = validate(search, reads, result);
            if (actual.digest != expected.digest || actual.pieces != expected.pieces)
                throw std::runtime_error("seed output changed between repeated sweeps");
            rounds.push_back({total_seconds/static_cast<double>(sweeps), total_seconds, cpu_seconds, actual.digest});
        }
        const auto peak_rss = peakRssBytes();
        std::cout << std::setprecision(12)
            << "{\n  \"schema\": \"nustar-seed-benchmark-v1\",\n"
            << "  \"reads\": " << reads.size() << ",\n"
            << "  \"read_limit\": " << limit << ",\n"
            << "  \"input_bases\": " << input_bases << ",\n"
            << "  \"paired\": " << (paired ? "true" : "false") << ",\n"
            << "  \"pieces_per_sweep\": " << expected.pieces << ",\n"
            << "  \"timed_sweeps\": " << sweeps << ",\n"
            << "  \"digest_algorithm\": \"sha256-full-ordered-result-v1\",\n"
            << "  \"digest\": \"" << expected.digest << "\",\n"
            << "  \"timing\": \"bulk reused-output searches; validation in separate untimed sweeps\",\n"
            << "  \"peak_rss_bytes\": " << peak_rss << ",\n"
            << "  \"policy\": {\"seedSearchLmax\": " << config.seedSearchLmax
            << ", \"seedSearchStartLmax\": " << config.seedSearchStartLmax
            << ", \"seedMultimapNmax\": " << config.seedMultimapNmax
            << ", \"seedSplitMin\": " << config.seedSplitMin
            << ", \"seedMapMin\": " << config.seedMapMin << "},\n"
            << "  \"rounds\": [\n";
        for (std::size_t i = 0; i < rounds.size(); ++i)
            std::cout << "    {\"seconds\": " << rounds[i].seconds
                << ", \"total_seconds\": " << rounds[i].total_seconds
                << ", \"cpu_seconds\": " << rounds[i].cpu_seconds
                << ", \"timed_sweeps\": " << sweeps
                << ", \"reads\": " << reads.size() << ", \"digest\": \""
                << rounds[i].digest << "\"}" << (i+1 == rounds.size() ? "\n" : ",\n");
        std::cout << "  ]\n}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "seed-benchmark: " << error.what() << '\n';
        return 1;
    }
}
