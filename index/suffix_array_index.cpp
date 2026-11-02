#include "suffix_array_index.hpp"
#include "sais.hpp"
#include "wide_suffix_array.hpp"
#include "../reference/genome.hpp"
#include "../common/binary.hpp"
#include "../common/file_identity.hpp"
#include "../common/sha256.hpp"
#include "star_rank.hpp"
#include "physical_sequence.hpp"

#include <algorithm>
#include <atomic>
#include <bit>
#include <charconv>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <exception>
#include <fstream>
#include <limits>
#include <stdexcept>

#if defined(__unix__)
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

namespace nustar::star2711b {
namespace {

// The immutable index artifact owns its mapping through every SA/SAi lease.
// Mapping avoids a second 26 GB allocation and copy. Like STAR, loading
// trusts index payloads after reading their format and size metadata.
struct MappedIndexFile {
    const std::uint8_t* data=nullptr;
    std::size_t size=0;
    ~MappedIndexFile() {
#if defined(__unix__)
        if (data) ::munmap(const_cast<std::uint8_t*>(data),size);
#endif
    }
    static std::shared_ptr<MappedIndexFile> open(const std::filesystem::path& path,
        std::uint64_t expected_size) {
#if defined(__unix__)
        auto result=std::make_shared<MappedIndexFile>();
        const int fd=::open(path.c_str(),O_RDONLY|O_CLOEXEC);
        if (fd<0) return {};
        struct stat st{};
        if (::fstat(fd,&st)!=0 || st.st_size<0 || std::uint64_t(st.st_size)!=expected_size) {
            ::close(fd); return {};
        }
        int flags=MAP_PRIVATE;
#if defined(MAP_POPULATE)
        // Load pages here, as STAR does when reading SA/SAi, rather than move
        // all first-touch I/O into the measured alignment stage.
        flags |= MAP_POPULATE;
#endif
        void* address=::mmap(nullptr,expected_size,PROT_READ,flags,fd,0);
        ::close(fd);
        if (address==MAP_FAILED) return {};
        result->data=static_cast<const std::uint8_t*>(address); result->size=expected_size;
        // Keep loaded pages available for random SA/SAi access during mapping.
        ::madvise(address,expected_size,MADV_NORMAL);
        return result;
#else
        return {};
#endif
    }
    void ready(std::uint64_t retained_bytes) const {
#if defined(__unix__)
        ::madvise(const_cast<std::uint8_t*>(data),size,MADV_RANDOM);
        // A v4 file without a trailing guard uses an owned SAi buffer. Release
        // its copied file-backed pages, retaining the SA's guard bytes.
        const auto page=std::uint64_t(::sysconf(_SC_PAGESIZE));
        const auto first=(retained_bytes+page-1)/page*page;
        if (first<size) ::madvise(const_cast<std::uint8_t*>(data)+first,size-first,MADV_DONTNEED);
#endif
    }
};

std::atomic<std::uint64_t> next_index_generation{1};

Failure indexError(ErrorKind kind, std::string message) {
    return Failure{kind, Stage::Reference, std::nullopt, 1, std::move(message)};
}

struct BuildParameters {
    std::uint32_t prefix_bases = 10;
    int runThreadN = 0;
};
constexpr std::uint32_t star_suffix_length_max =
    std::numeric_limits<std::uint32_t>::max();

Result<BuildParameters> buildParameters(const BackendOptions& options) {
    if (!options.schema_id.empty() &&
        !((options.schema_id == "nustar-sa-v2" && options.schema_version == 2) ||
          (options.schema_id == "nustar-sa-v3" && options.schema_version == 3) ||
          (options.schema_id == "nustar-sa-v4" && options.schema_version == 4)))
        return indexError(ErrorKind::Config, "unsupported suffix array options schema");
    BuildParameters parameters;
    for (const auto& option : options.values) {
        if (option.resolved_tokens.size() != 1)
            return indexError(ErrorKind::Config, "invalid suffix array option: " + option.name);
        const auto& token = option.resolved_tokens[0];
        std::uint32_t value = 0;
        const auto [end, error] = std::from_chars(token.data(), token.data()+token.size(), value);
        if (error != std::errc{} || end != token.data()+token.size() || value == 0)
            return indexError(ErrorKind::Config, "invalid suffix array option: " + option.name);
        if (option.name == "genomeSAindexNbases" && value <= 15)
            parameters.prefix_bases = value;
        else if (option.name == "runThreadN" && value <= 256)
            parameters.runThreadN = static_cast<int>(value);
        else return indexError(ErrorKind::Config, "unknown or out-of-range suffix array option: " +
                               option.name);
    }
    return parameters;
}

std::string optionsHash(const BuildParameters& parameters, std::string_view algorithm,
                        std::uint32_t format_version) {
    Sha256 hash;
    hash.update("nustar-sa-v"+std::to_string(format_version)+":");
    hash.update(std::to_string(parameters.prefix_bases));
    hash.update(":");
    hash.update(std::to_string(star_suffix_length_max));
    hash.update(":");
    hash.update(algorithm);
    return hash.finish();
}

template<class Position>
std::vector<PrefixRange> buildPrefixes(std::span<const std::uint8_t> sequence,
    std::span<const Position> suffixes, std::uint32_t bases, int runThreadN) {
    const auto domain = std::uint64_t{1} << (2*bases);
    auto scan = [&](std::uint64_t first, std::uint64_t last, std::uint64_t reserve) {
        std::vector<PrefixRange> result;
        result.reserve(std::min(reserve, last-first));
        for (auto i = first; i < last; ++i) {
            const auto position = suffixes[i];
            if (position + bases > sequence.size()) continue;
            std::uint64_t code = 0;
            bool valid = true;
            for (std::uint32_t j = 0; j < bases; ++j) {
                const auto symbol = sequence[position+j];
                if (symbol < 1 || symbol > 4) { valid = false; break; }
                code = (code << 2) | (symbol - 1);
            }
            if (!valid) continue;
            if (result.empty() || result.back().code != code)
                result.push_back(PrefixRange{code, i, i+1});
            else result.back().end = i+1;
        }
        return result;
    };
    int threads = 1;
#if defined(_OPENMP)
    threads = detail::suffixConstructionThreads(sequence.size(), runThreadN);
#endif
    if (threads == 1) return scan(0, suffixes.size(), domain);
    // SA chunks preserve rank order. Only a prefix spanning a chunk boundary
    // needs merging; every entry and range is identical to the serial scan.
    std::vector<std::vector<PrefixRange>> chunks(threads);
    std::vector<std::exception_ptr> errors(threads);
#if defined(_OPENMP)
    #pragma omp parallel for schedule(static) num_threads(threads)
#endif
    for (int part = 0; part < threads; ++part) {
        const auto size = suffixes.size()/threads;
        const auto extra = suffixes.size()%threads;
        const auto first = size*part+std::min<std::size_t>(part, extra);
        const auto last = size*(part+1)+std::min<std::size_t>(part+1, extra);
        try { chunks[part] = scan(first, last, domain/threads+1); }
        catch (...) { errors[part] = std::current_exception(); }
    }
    std::size_t count = 0;
    for (int part = 0; part < threads; ++part) {
        if (errors[part]) std::rethrow_exception(errors[part]);
        count += chunks[part].size();
    }
    std::vector<PrefixRange> result;
    result.reserve(count);
    for (auto& chunk : chunks) {
        auto first = chunk.begin();
        if (!result.empty() && first != chunk.end() && result.back().code == first->code) {
            result.back().end = first->end;
            ++first;
        }
        result.insert(result.end(), first, chunk.end());
        std::vector<PrefixRange>().swap(chunk);
    }
    return result;
}

void writeIndex(const std::filesystem::path& path, const SuffixArrayIndex& index) {
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) throw std::runtime_error("cannot create index: " + path.string());
    out.write("KSAI", 4);
    binary::put32(out, 4);
    binary::put32(out, index.prefixBases());
    binary::put32(out, index.suffixLengthMax());
    binary::put64(out, index.physicalLength());
    binary::put64(out, index.suffixArray().size());
    binary::put64(out, StarPrefixIndex::entryCount(index.prefixBases()));
    binary::putString(out, index.identity().reference_semantic_id);
    binary::putString(out, optionsHash(BuildParameters{index.prefixBases()}, index.algorithmId(), 4));
    binary::putString(out, std::string(index.algorithmId()));
    binary::put32(out, index.suffixArray().bitWidth());
    const auto packed = index.suffixArray().packedBytes();
    out.write(reinterpret_cast<const char*>(packed.data()),
              static_cast<std::streamsize>(packed.size()));
    const auto prefix = index.prefixIndex().bytes();
    out.write(reinterpret_cast<const char*>(prefix.data()),
              static_cast<std::streamsize>(prefix.size()));
    const char guard[7]{};
    out.write(guard,sizeof(guard));
    out.close();
    if (!out) throw std::runtime_error("failed writing index: " + path.string());
}

void writeChecksum(const std::filesystem::path& path, std::string_view checksum) {
    std::ofstream out(path.string() + ".sha256", std::ios::trunc);
    if (!out || !(out << checksum << '\n'))
        throw std::runtime_error("cannot write index checksum: " + path.string());
}

} // namespace

SuffixArrayIndex::SuffixArrayIndex(SearchIndexIdentity identity, ReferenceLease reference,
    std::filesystem::path artifact, std::string algorithm_id,
    std::uint32_t prefix_bases, std::uint32_t suffix_length_max,
    std::uint64_t physical_length, PackedSuffixArray suffixes,
    std::vector<PrefixRange> prefixes, StarPrefixIndex native_prefix)
    : identity_(std::move(identity)), reference_(std::move(reference)),
      artifact_(std::move(artifact)), algorithm_id_(std::move(algorithm_id)),
      prefix_bases_(prefix_bases), suffix_length_max_(suffix_length_max),
      physical_length_(physical_length), suffixes_(std::move(suffixes)),
      prefix_index_(std::move(native_prefix)) {
    if (prefix_index_.bytes().empty()) {
        const auto& genome=dynamic_cast<const GenomeReferenceData&>(*reference_);
        prefix_index_=StarPrefixIndex(prefix_bases_,
            genome.storage().forward.size()+genome.synthetic().size(),
            suffixes_,genome,prefixes);
    }
}

SearchIndexCapabilities SuffixArrayIndex::capabilities() const {
    return SearchIndexCapabilities{{"nustar-sa-exact-v2"}, false};
}

SuffixArrayIndexBackend::SuffixArrayIndexBackend(SuffixArrayBuilderLease builder)
    : builder_(builder ? std::move(builder) : std::make_shared<SaisSuffixArrayBuilder>()) {
    if (builder_->algorithmId().empty())
        throw std::invalid_argument("suffix array builder needs an algorithm ID");
}

Result<SearchIndexLease> SuffixArrayIndexBackend::build(const SearchIndexRequest& request,
    ReferenceLease reference) {
    if (request.backend_id != backendId() || request.artifact.empty())
        return indexError(ErrorKind::Config, "suffix array backend or artifact path is invalid");
    const auto concrete = std::dynamic_pointer_cast<const GenomeReferenceData>(reference);
    if (!concrete) return indexError(ErrorKind::Config, "suffix array requires a NuSTAR reference");
    auto option_result = buildParameters(request.options);
    if (std::holds_alternative<Failure>(option_result))
        return std::get<Failure>(std::move(option_result));
    const auto options = std::get<BuildParameters>(option_result);
    std::vector<std::uint8_t> sequence;
    PackedSuffixArray packed;
    std::vector<PrefixRange> prefixes;
    try {
        sequence = makeStarPhysicalSequence(*concrete, 1);
        auto finish = [&]<class Position>(std::vector<Position> all) {
            if (all.size() != sequence.size()+1)
                throw std::runtime_error("suffix array builder returned the wrong length");
            std::size_t kept = 0;
            for (const auto position : all)
                if (position < sequence.size() && sequence[position] >= 1 && sequence[position] <= 4)
                    all[kept++] = position;
            all.resize(kept);
            if constexpr (sizeof(Position) == 8)
                detail::normalizeStarSuffixRanks64(all, sequence, concrete->storage().forward.size(),
                    concrete->synthetic().size(), star_suffix_length_max, options.runThreadN);
            else
                normalizeStarSuffixRanks(all, sequence, concrete->storage().forward.size(),
                    concrete->synthetic().size(), star_suffix_length_max);
            prefixes = buildPrefixes(sequence, std::span<const Position>(all), options.prefix_bases, options.runThreadN);
            packed = PackedSuffixArray(all, sequence.size());
        };
        if (dynamic_cast<const SaisSuffixArrayBuilder*>(builder_.get()))
            finish(detail::buildSuffixArraySais64(sequence, options.runThreadN));
        else {
            // Keep existing caller-supplied 32-bit builder implementations ABI/source-compatible.
            if (sequence.size() >= std::numeric_limits<std::uint32_t>::max())
                throw std::length_error("configured custom builder only supports 32-bit positions");
            finish(builder_->build(sequence));
        }
    } catch (const std::exception& error) {
        return indexError(ErrorKind::Capacity, error.what());
    }
    auto identity = SearchIndexIdentity{"nustar-sa", "nustar-sa-v4", 4, "",
        std::string(concrete->semanticId()), optionsHash(options, builderId(), 4),
        IndexGeneration{next_index_generation.fetch_add(1)}};
    std::error_code ec;
    if (!request.artifact.parent_path().empty())
        std::filesystem::create_directories(request.artifact.parent_path(), ec);
    if (ec) return indexError(ErrorKind::Io, "cannot create index directory: " + ec.message());
    auto index = std::make_shared<SuffixArrayIndex>(identity, reference, request.artifact,
        std::string(builderId()), options.prefix_bases, star_suffix_length_max,
        sequence.size(), std::move(packed), std::move(prefixes));
    try { writeIndex(request.artifact, *index); }
    catch (const std::exception& error) { return indexError(ErrorKind::Io, error.what()); }
    auto file = identifyFile(request.artifact);
    if (std::holds_alternative<Failure>(file)) return std::get<Failure>(std::move(file));
    index->setArtifactSha256(std::get<FileIdentity>(file).sha256);
    try { writeChecksum(request.artifact, index->identity().artifact_sha256); }
    catch (const std::exception& error) { return indexError(ErrorKind::Io, error.what()); }
    return SearchIndexLease{std::move(index)};
}

Result<SearchIndexLease> SuffixArrayIndexBackend::load(const SearchIndexRequest& request,
    ReferenceLease reference) {
    if (request.backend_id != backendId())
        return indexError(ErrorKind::Config, "wrong index backend");
    const auto concrete = std::dynamic_pointer_cast<const GenomeReferenceData>(reference);
    if (!concrete) return indexError(ErrorKind::Config, "suffix array requires a NuSTAR reference");
    std::ifstream in(request.artifact, std::ios::binary);
    if (!in) return indexError(ErrorKind::Io, "cannot open index: " + request.artifact.string());
    char magic[4]{}; in.read(magic, 4);
    std::uint32_t version=0, bases=0, suffix_length_max=0, packed_width=0;
    std::uint64_t physical=0, suffix_count=0, prefix_count=0;
    std::string semantic, options_hash, algorithm_id;
    if (std::string_view(magic,4)!="KSAI" || !binary::get32(in,version) ||
        (version!=2 && version!=3 && version!=4) ||
        !binary::get32(in,bases) || bases<1 || bases>15 ||
        !binary::get32(in,suffix_length_max) ||
            suffix_length_max != star_suffix_length_max ||
        !binary::get64(in,physical) || !binary::get64(in,suffix_count) ||
        !binary::get64(in,prefix_count) ||
        !binary::getString(in,semantic,128) || !binary::getString(in,options_hash,128) ||
        !binary::getString(in,algorithm_id,128) || algorithm_id.empty() ||
        semantic!=concrete->semanticId())
        return indexError(ErrorKind::Index, "invalid index header or reference identity");
    if (version >= 3 && !binary::get32(in, packed_width))
        return indexError(ErrorKind::Index, "truncated packed suffix array header");
    // The stored option digest is provenance, not a format discriminator.
    // Dimensions and suffix_length_max above describe the actual encoding;
    // changing the product name must not invalidate an otherwise identical SA.
    if (!request.options.schema_id.empty()) {
        auto requested = buildParameters(request.options);
        if (std::holds_alternative<Failure>(requested))
            return std::get<Failure>(std::move(requested));
        const auto expected = std::get<BuildParameters>(requested);
        if (bases != expected.prefix_bases)
            return indexError(ErrorKind::Index, "index build options do not match request");
    }
    std::error_code ec;
    const auto bytes = std::filesystem::file_size(request.artifact, ec);
    const auto header_bytes = in.tellg();
    const auto forward_size = concrete->storage().forward.size();
    const auto synthetic_size = concrete->synthetic().size();
    constexpr auto max_reference_bases =
        (std::numeric_limits<std::uint64_t>::max()-80)/2;
    if (forward_size > max_reference_bases ||
        synthetic_size > max_reference_bases-forward_size)
        return indexError(ErrorKind::Capacity, "reference length exceeds index format");
    const std::uint64_t expected_physical =
        2 * (std::uint64_t(forward_size) + std::uint64_t(synthetic_size)) + 80;
    if (ec || header_bytes < 0 || physical != expected_physical ||
        physical == 0 ||
        (version == 2 && physical >= static_cast<std::uint64_t>(
            std::numeric_limits<std::int32_t>::max())) ||
        suffix_count > physical || (version<4 && prefix_count > suffix_count) ||
        (version==4 && prefix_count!=StarPrefixIndex::entryCount(bases)) ||
        (version >= 3 && packed_width != PackedSuffixArray::widthFor(physical)))
        return indexError(ErrorKind::Index, "invalid index lengths or packed width");
    if (suffix_count > std::numeric_limits<std::size_t>::max() ||
        prefix_count > std::numeric_limits<std::size_t>::max() /
            sizeof(PrefixRange))
        return indexError(ErrorKind::Index, "index lengths exceed addressable memory");
    std::uint64_t suffix_bytes = 0;
    if (version == 2) {
        if (suffix_count > std::numeric_limits<std::uint64_t>::max()/4)
            return indexError(ErrorKind::Index, "legacy suffix array size overflow");
        suffix_bytes = 4*suffix_count;
    } else {
        try {
            suffix_bytes = PackedSuffixArray::requiredBytes(
                static_cast<std::size_t>(suffix_count),
                static_cast<std::uint8_t>(packed_width));
        } catch (const std::invalid_argument&) {
            return indexError(ErrorKind::Index, "packed suffix array size overflow");
        }
    }
    const auto prefix_bytes=version==4 ? PackedSuffixArray::requiredBytes(
        prefix_count,StarPrefixIndex::width(forward_size+synthetic_size)) : 24*prefix_count;
    const auto expected_bytes=static_cast<std::uint64_t>(header_bytes)+suffix_bytes+prefix_bytes;
    const bool prefix_guard=version==4 && bytes==expected_bytes+7;
    if (header_bytes < 0 ||
        static_cast<std::uint64_t>(header_bytes) > bytes ||
        suffix_bytes > bytes-static_cast<std::uint64_t>(header_bytes) ||
        prefix_bytes > bytes-static_cast<std::uint64_t>(header_bytes)-suffix_bytes ||
        (bytes != expected_bytes && !prefix_guard))
        return indexError(ErrorKind::Index, "invalid index lengths");
    if (suffix_bytes > static_cast<std::uint64_t>(
            std::numeric_limits<std::streamsize>::max()) ||
        prefix_count > static_cast<std::uint64_t>(
            std::numeric_limits<std::streamsize>::max()) / sizeof(PrefixRange))
        return indexError(ErrorKind::Index, "index payload exceeds stream limits");
    using Clock=std::chrono::steady_clock;
    const auto read_started=Clock::now();
    const auto mapping=version==4 ? MappedIndexFile::open(request.artifact,bytes) : nullptr;
    PackedSuffixArray suffixes;
    try {
        if (mapping) {
            suffixes=PackedSuffixArray::fromView(packed_width,suffix_count,
                std::span(mapping->data,mapping->size).subspan(static_cast<std::size_t>(header_bytes)),mapping,true);
            in.seekg(suffix_bytes,std::ios::cur);
            if (!in) return indexError(ErrorKind::Index,"invalid mapped SA offset");
        } else if (version >= 3) {
            // fromBytes appends seven bytes for safe unaligned word loads.
            // Reserve them before reading: growing an exact-size human SA
            // would allocate and copy another ~26 GB merely to add padding.
            std::vector<std::uint8_t> payload(static_cast<std::size_t>(suffix_bytes)+7);
            payload.resize(static_cast<std::size_t>(suffix_bytes));
            in.read(reinterpret_cast<char*>(payload.data()),
                    static_cast<std::streamsize>(payload.size()));
            if (!in) return indexError(ErrorKind::Index, "truncated packed suffix array");
            suffixes = PackedSuffixArray::fromBytes(
                static_cast<std::uint8_t>(packed_width),
                static_cast<std::size_t>(suffix_count), std::move(payload));
        } else {
            std::vector<std::uint32_t> legacy(static_cast<std::size_t>(suffix_count));
            if constexpr (std::endian::native == std::endian::little)
                in.read(reinterpret_cast<char*>(legacy.data()),
                        static_cast<std::streamsize>(legacy.size() * sizeof(std::uint32_t)));
            else for (auto& position : legacy)
                if (!binary::get32(in, position))
                    return indexError(ErrorKind::Index, "truncated legacy suffix array");
            if (!in) return indexError(ErrorKind::Index, "truncated legacy suffix array");
            suffixes = PackedSuffixArray(legacy, physical);
        }
    } catch (const std::exception& error) {
        return indexError(ErrorKind::Index, error.what());
    }
    const auto prefix_started=Clock::now();
    std::vector<PrefixRange> prefixes;
    StarPrefixIndex native_prefix;
    if (version==4) {
        try {
            if (mapping && prefix_guard) {
                native_prefix=StarPrefixIndex::fromView(bases,forward_size+synthetic_size,
                    std::span(mapping->data,mapping->size).subspan(
                        static_cast<std::size_t>(header_bytes)+suffix_bytes),mapping);
            } else {
                std::vector<std::uint8_t> payload(prefix_bytes+7);
                payload.resize(prefix_bytes);
                in.read(reinterpret_cast<char*>(payload.data()),payload.size());
                if (!in) return indexError(ErrorKind::Index,"truncated STAR prefix index");
                native_prefix=StarPrefixIndex::fromBytes(bases,forward_size+synthetic_size,
                    std::move(payload));

            }
        } catch (const std::exception& e) { return indexError(ErrorKind::Index,e.what()); }
    } else {
    prefixes.resize(static_cast<std::size_t>(prefix_count));
    if constexpr (std::endian::native == std::endian::little)
        in.read(reinterpret_cast<char*>(prefixes.data()),
                static_cast<std::streamsize>(prefixes.size() * sizeof(PrefixRange)));
    else for (auto& entry : prefixes)
        if (!binary::get64(in,entry.code) || !binary::get64(in,entry.begin) ||
            !binary::get64(in,entry.end))
            return indexError(ErrorKind::Index, "truncated prefix index");
    if (!in) return indexError(ErrorKind::Index, "truncated prefix index");
    }
    // Optional build provenance only; STAR does not rehash index payloads on
    // load. A missing or stale sidecar never prevents ordinary alignment.
    std::ifstream provenance(request.artifact.string() + ".sha256");
    std::string artifact_digest;
    std::getline(provenance, artifact_digest);
    if (std::getenv("NUSTAR_INDEX_LOAD_TIMING")) {
        const auto seconds=[](auto a, auto b) { return std::chrono::duration<double>(b-a).count(); };
        std::cerr << "index_load_seconds sa_read=" << seconds(read_started,prefix_started)
            << " prefix_read=" << seconds(prefix_started,Clock::now()) << '\n';
    }
    if (mapping) mapping->ready(prefix_guard ? bytes :
        static_cast<std::uint64_t>(header_bytes)+suffix_bytes+7);
    auto identity = SearchIndexIdentity{"nustar-sa",
        "nustar-sa-v"+std::to_string(version), version,
        std::move(artifact_digest), semantic, options_hash,
        IndexGeneration{next_index_generation.fetch_add(1)}};
    return SearchIndexLease{std::make_shared<const SuffixArrayIndex>(std::move(identity),
        std::move(reference), request.artifact, algorithm_id, bases, suffix_length_max, physical,
        std::move(suffixes), std::move(prefixes), std::move(native_prefix))};
}

Status SuffixArrayIndexBackend::save(const ISearchIndex& source,
    const std::filesystem::path& target) {
    const auto* index = dynamic_cast<const SuffixArrayIndex*>(&source);
    if (!index) return indexError(ErrorKind::Config, "cannot save another index backend");
    // Export the in-memory layout as v4, including indexes loaded from v2/v3.
    // A separate temporary file ensures failed writes never replace the source.
    const auto temporary=std::filesystem::path(target.string()+".writing");
    try {
        writeIndex(temporary,*index);
        auto file=identifyFile(temporary);
        if (std::holds_alternative<Failure>(file)) return std::get<Failure>(std::move(file));
        const auto hash=std::get<FileIdentity>(file).sha256;
        std::filesystem::rename(temporary,target);
        writeChecksum(target,hash);
    } catch (const std::exception& error) { return indexError(ErrorKind::Io,error.what()); }
    return std::monostate{};
}

Result<SearchIndexLease> SuffixArrayIndexBackend::updateAnnotations(
    SearchIndexLease, ReferenceLease, const BackendOptions&) {
    return indexError(ErrorKind::Config,
        "suffix array backend requires a full rebuild for annotation updates");
}

} // namespace nustar::star2711b
