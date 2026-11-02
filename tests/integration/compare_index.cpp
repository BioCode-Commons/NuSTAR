#include "../../reference/genome.hpp"
#include "../../common/binary.hpp"
#include "../../index/suffix_array_index.hpp"

#include <algorithm>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <vector>

using namespace nustar::star2711b;

namespace {

template<class T> T checked(Result<T> result) {
    if (std::holds_alternative<Failure>(result))
        throw std::runtime_error(std::get<Failure>(result).message);
    return std::get<T>(std::move(result));
}

std::vector<std::uint8_t> readBytes(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary | std::ios::ate);
    if (!input) throw std::runtime_error("cannot open " + path.string());
    const auto length = input.tellg();
    if (length < 0) throw std::runtime_error("cannot size " + path.string());
    std::vector<std::uint8_t> bytes(static_cast<std::size_t>(length));
    input.seekg(0);
    if (!input.read(reinterpret_cast<char*>(bytes.data()), length))
        throw std::runtime_error("cannot read " + path.string());
    return bytes;
}

std::uint64_t packed(std::span<const std::uint8_t> bytes, std::uint64_t rank,
                     unsigned word_bits) {
    const auto bit = rank * word_bits;
    const auto offset = bit / 8;
    std::uint64_t word = 0;
    for (unsigned i = 0; i < 8; ++i) word |= std::uint64_t(bytes[offset+i]) << (8*i);
    return (word >> (bit % 8)) & ((std::uint64_t{1} << word_bits) - 1);
}

std::uint64_t starValue(std::uint64_t position, std::uint64_t real,
                        std::uint64_t synthetic, unsigned strand_bit) {
    const auto forward_synthetic = 2*real+40;
    const auto reverse_synthetic = forward_synthetic+synthetic+20;
    if (position < real) return position;
    if (position >= real+20 && position < 2*real+20)
        return (std::uint64_t{1} << strand_bit) |
               (synthetic + position - (real+20));
    if (position >= forward_synthetic && position < forward_synthetic+synthetic)
        return real + position - forward_synthetic;
    if (position >= reverse_synthetic && position < reverse_synthetic+synthetic)
        return (std::uint64_t{1} << strand_bit) | (position - reverse_synthetic);
    throw std::runtime_error("NuSTAR SA contains a separator position");
}

std::optional<std::uint64_t> starPrefix(std::uint64_t value, unsigned strand_bit,
    std::span<const std::uint8_t> genome, unsigned bases) {
    const auto strand = (value >> strand_bit) != 0;
    const auto offset = value & ((std::uint64_t{1} << strand_bit) - 1);
    std::uint64_t code = 0;
    for (unsigned i = 0; i < bases; ++i) {
        if (offset+i >= genome.size()) return {};
        auto base = genome[strand ? genome.size()-1-offset-i : offset+i];
        if (base > 3) return {};
        if (strand) base = 3-base;
        code = (code << 2) | base;
    }
    return code;
}

std::uint8_t starBase(std::uint64_t value, unsigned strand_bit,
    std::span<const std::uint8_t> genome, std::uint64_t distance) {
    const auto strand = (value >> strand_bit) != 0;
    const auto offset = value & ((std::uint64_t{1} << strand_bit) - 1);
    if (offset+distance >= genome.size()) return 6;
    auto base = genome[strand ? genome.size()-1-offset-distance : offset+distance];
    if (strand && base < 4) base = 3-base;
    return base;
}

std::optional<std::uint64_t> nustarPrefix(std::uint64_t position,
    std::span<const std::uint8_t> genome, std::uint64_t real,
    std::uint64_t synthetic, unsigned bases) {
    const auto forward_synthetic = 2*real+40;
    const auto reverse_synthetic = forward_synthetic+synthetic+20;
    std::uint64_t code = 0;
    for (unsigned i = 0; i < bases; ++i) {
        const auto p = std::uint64_t(position)+i;
        std::uint8_t base = 5;
        if (p < real) base = genome[p];
        else if (p >= real+20 && p < 2*real+20) {
            base = genome[real-1-(p-(real+20))];
            if (base < 4) base = 3-base;
        } else if (p >= forward_synthetic && p < forward_synthetic+synthetic)
            base = genome[real+p-forward_synthetic];
        else if (p >= reverse_synthetic && p < reverse_synthetic+synthetic) {
            base = genome[genome.size()-1-(p-reverse_synthetic)];
            if (base < 4) base = 3-base;
        }
        if (base > 3) return {};
        code = (code << 2) | base;
    }
    return code;
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 3 || argc > 4 ||
        (argc == 4 && std::string_view(argv[3]) != "--show-rank-diffs")) {
        std::cerr << "Usage: compare-genome-indices STAR_DIR NUSTAR_DIR [--show-rank-diffs]\n";
        return 2;
    }
    try {
        const std::filesystem::path star_dir(argv[1]), nustar_dir(argv[2]);
        auto reference = checked(loadGenomeReference(nustar_dir));
        SearchIndexRequest request{"nustar-sa", nustar_dir / "index.sa", {}};
        SuffixArrayIndexBackend backend;
        auto index = std::dynamic_pointer_cast<const SuffixArrayIndex>(
            checked(backend.load(request, reference)));
        if (!index) throw std::runtime_error("loaded index is not a suffix array");
        const auto real = reference->storage().forward.size();
        const auto synthetic = reference->synthetic().size();
        const auto genome = readBytes(star_dir / "Genome");
        if (genome.size() != real+synthetic)
            throw std::runtime_error("STAR Genome length differs");
        for (std::size_t i = 0; i < real; ++i)
            if (genome[i] != static_cast<std::uint8_t>(reference->storage().forward[i]))
                throw std::runtime_error("STAR Genome forward sequence differs");
        for (std::size_t i = 0; i < synthetic; ++i)
            if (genome[real+i] != static_cast<std::uint8_t>(reference->synthetic()[i]))
                throw std::runtime_error("STAR Genome synthetic sequence differs");
        const auto sa = readBytes(star_dir / "SA");
        const auto& suffixes = index->suffixArray();
        if (suffixes.empty()) throw std::runtime_error("empty NuSTAR suffix array");
        unsigned strand_bit = 0;
        for (unsigned bit = 32; bit < 63; ++bit)
            if ((suffixes.size()-1)*(bit+1)/8+8 == sa.size()) {
                strand_bit = bit;
                break;
            }
        if (!strand_bit || genome.size() >= (std::uint64_t{1} << strand_bit))
            throw std::runtime_error("STAR packed SA size or strand bit differs");

        // Compare every native SAi entry, including absent and N flags at
        // every prefix length, independently of final mapping output.
        const auto native_bytes=readBytes(star_dir / "SAindex");
        std::ifstream native_file(star_dir / "SAindex",std::ios::binary);
        std::uint64_t native_bases=0;
        if (!binary::get64(native_file,native_bases) || native_bases!=index->prefixBases())
            throw std::runtime_error("STAR prefix depth differs");
        const auto& prefix=index->prefixIndex();
        for (unsigned level=0;level<=native_bases;++level) {
            std::uint64_t offset=0;
            if (!binary::get64(native_file,offset) ||
                offset!=(level ? prefix.levelEnd(level) : 0))
                throw std::runtime_error("STAR prefix level offset differs");
        }
        const auto payload=std::span(native_bytes).subspan(8*(native_bases+2));
        for (std::uint64_t i=0;i<StarPrefixIndex::entryCount(native_bases);++i)
            {
                const auto star_entry=packed(payload,i,strand_bit+3), ours=prefix[i];
                // sjdbBuildIndex leaves obsolete ranks in absent entries at
                // the end of a level. Search tests absent before reading rank.
                const bool absent=star_entry&prefix.absentMask();
                if (bool(ours&prefix.absentMask())!=absent ||
                    (!absent && star_entry!=ours))
                    throw std::runtime_error("STAR SAi entry differs at "+std::to_string(i));
            }
        std::cout << "SAi entries: " << StarPrefixIndex::entryCount(native_bases)
                  << ", all presence flags and present ranks/N flags identical\n";

        // Exact equality at every rank already proves membership/permutation
        // equality to STAR. No random-access seen[] pass or sampled ranks.
        const auto ranges = index->prefixRanges();
        std::uint64_t previous_end = 0, previous_code = 0;
        bool have_previous = false;
        for (const auto& range : ranges) {
            if (range.begin >= range.end || range.end > suffixes.size() ||
                range.begin < previous_end || (have_previous && range.code <= previous_code))
                throw std::runtime_error("invalid or overlapping prefix ranges");
            previous_end = range.end;
            previous_code = range.code;
            have_previous = true;
        }
        struct alignas(64) Check {
            std::uint64_t exact=0, prefix=0, ranges=0, maximum=0, above=0;
            std::uint64_t padding=0, separator=0, unexpected=0;
            std::vector<std::uint64_t> shorter;
            std::vector<std::string> diagnostics;
            std::exception_ptr error;
        };
        int threads = 1;
#if defined(_OPENMP)
        if (suffixes.size() >= (1ULL << 20)) threads = 16;
#endif
        std::vector<Check> checks(threads);
#if defined(_OPENMP)
        #pragma omp parallel for schedule(static) num_threads(threads)
#endif
        for (int part = 0; part < threads; ++part) {
          auto& check = checks[part];
          try {
            check.shorter.resize(index->prefixBases()+1);
            const auto size = suffixes.size()/threads, extra = suffixes.size()%threads;
            const auto first = size*part+std::min<std::size_t>(part, extra);
            const auto last = size*(part+1)+std::min<std::size_t>(part+1, extra);
            auto range = std::lower_bound(ranges.begin(), ranges.end(), first,
                [](const auto& entry, auto rank) { return entry.end <= rank; });
            for (auto rank = first; rank < last; ++rank) {
                const auto position = suffixes[rank];
                check.maximum = std::max(check.maximum, position);
                check.above += position > UINT32_MAX;
                const auto star_value = packed(sa, rank, strand_bit+1);
                const auto mapped = starValue(position, real, synthetic, strand_bit);
                if (mapped == star_value) ++check.exact;
                else {
                    std::uint64_t distance = 0;
                    for (; distance < index->suffixLengthMax(); ++distance) {
                        const auto a = starBase(mapped, strand_bit, genome, distance);
                        const auto b = starBase(star_value, strand_bit, genome, distance);
                        if (a != b || a >= 5) break;
                    }
                    const auto a = starBase(mapped, strand_bit, genome, distance);
                    const auto b = starBase(star_value, strand_bit, genome, distance);
                    if (distance < index->suffixLengthMax() && a == b && a >= 5)
                        ++check.padding;
                    else if ((a == 5 && b == 6) || (a == 6 && b == 5))
                        ++check.separator;
                    else ++check.unexpected;
                    if (argc == 4 && check.diagnostics.size() < 14) {
                        std::ostringstream message;
                        message << "rank " << rank << ": NuSTAR=" << mapped << ", STAR=" << star_value
                                << ", shared=" << distance << ", next=" << unsigned(a) << "/" << unsigned(b);
                        check.diagnostics.push_back(message.str());
                    }
                }
                const auto kp = nustarPrefix(position, genome, real, synthetic, index->prefixBases());
                const auto sp = starPrefix(star_value, strand_bit, genome, index->prefixBases());
                check.prefix += kp != sp;
                if (!kp && !sp)
                    for (unsigned length = 1; length < index->prefixBases(); ++length)
                        check.shorter[length] += nustarPrefix(position, genome, real, synthetic, length) !=
                                                 starPrefix(star_value, strand_bit, genome, length);
                while (range != ranges.end() && range->end <= rank) ++range;
                const bool covered = range != ranges.end() && range->begin <= rank;
                // Check every rank's membership, independently of the builder's
                // chunk/merge algorithm. Invalid prefixes must remain uncovered.
                check.ranges += kp ? (!covered || range->code != *kp) : covered;
            }
          } catch (...) { check.error = std::current_exception(); }
        }
        std::uint64_t maximum_position = 0, positions_above_32bit = 0;
        std::uint64_t exact_rank = 0, prefix_mismatch = 0;
        std::uint64_t padding_ties = 0, separator_differences = 0,
                      unexpected_rank_differences = 0;
        std::uint64_t range_mismatch = 0;
        unsigned shown_differences = 0;
        std::vector<std::uint64_t> shorter_mismatch(index->prefixBases()+1);
        for (const auto& check : checks) {
            if (check.error) std::rethrow_exception(check.error);
            maximum_position = std::max(maximum_position, check.maximum);
            positions_above_32bit += check.above;
            exact_rank += check.exact;
            prefix_mismatch += check.prefix;
            range_mismatch += check.ranges;
            padding_ties += check.padding;
            separator_differences += check.separator;
            unexpected_rank_differences += check.unexpected;
            for (unsigned length = 1; length < index->prefixBases(); ++length)
                shorter_mismatch[length] += check.shorter[length];
            for (const auto& message : check.diagnostics)
                if (shown_differences++ < 14) std::cerr << message << '\n';
        }
        std::cout << "SA entries: " << suffixes.size() << ", exact rank matches: "
                  << exact_rank << ", " << index->prefixBases()
                  << "-mer rank mismatches: " << prefix_mismatch
                  << ", padding ties: " << padding_ties
                  << ", separator differences: " << separator_differences
                  << ", unexpected rank differences: " << unexpected_rank_differences << '\n';
        std::cout << "Prefix range membership mismatches: " << range_mismatch << '\n';
        std::cout << "Position range: maximum=" << maximum_position
                  << ", above UINT32_MAX=" << positions_above_32bit
                  << ", packed bits=" << unsigned(suffixes.bitWidth())
                  << ", payload bytes=" << suffixes.packedBytes().size() << '\n';
        for (unsigned length = 1; length < index->prefixBases(); ++length)
            if (shorter_mismatch[length])
                std::cout << length << "-mer rank mismatches: "
                          << shorter_mismatch[length] << '\n';
        if (exact_rank != suffixes.size() || prefix_mismatch || range_mismatch || unexpected_rank_differences ||
            std::ranges::any_of(shorter_mismatch, [](auto count) { return count != 0; }))
            throw std::runtime_error("STAR and NuSTAR exact SA ranks differ");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "INDEX DIFFERENCE: " << error.what() << '\n';
        return 1;
    }
}
