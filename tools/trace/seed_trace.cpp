#include "../../seed_search/star_search.hpp"
#include "../../profiling/stage_progress.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace nustar::star2711b;

namespace {

struct FastqRecord {
    std::string name, sequence;
};

bool readFastq(std::ifstream& input, FastqRecord& record) {
    std::string plus, quality;
    if (!std::getline(input, record.name)) return false;
    if (!std::getline(input, record.sequence) || !std::getline(input, plus) ||
        !std::getline(input, quality) || record.name.empty() ||
        record.name[0] != '@' || plus.empty() || plus[0] != '+' ||
        quality.size() != record.sequence.size())
        throw std::runtime_error("invalid FASTQ record");
    const auto name_end = record.name.find_first_of(" \t");
    if (name_end != std::string::npos) record.name.resize(name_end);
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

MappingRead prepare(const FastqRecord& first, const FastqRecord* second,
                    std::uint64_t ordinal) {
    MappingRead read{};
    read.fragment = FragmentId{0, ordinal};
    auto& sequence = read.encoded[0];
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
    for (auto base : sequence) read.encoded[1].push_back(complement(base));
    read.encoded[2].assign(read.encoded[1].rbegin(), read.encoded[1].rend());
    read.budget = ReadBudget{0, sequence.size(), sequence.size(),
                            first.sequence.size()+(second ? second->sequence.size() : 0)};
    return read;
}

template<class T> T checked(Result<T> result) {
    if (std::holds_alternative<Failure>(result))
        throw std::runtime_error(std::get<Failure>(result).message);
    return std::get<T>(std::move(result));
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: seed-trace INDEX_DIR FIRST.fastq [SECOND.fastq] OUT.tsv "
                     "[--seedSearchLmax N] [--seedSearchStartLmax N] "
                     "[--seedMultimapNmax N] [--seedSplitMin N] [--seedMapMin N]\n";
        return 2;
    }
    try {
        const std::filesystem::path directory = argv[1];
        const bool paired = argc >= 5 && argv[4][0] != '-';
        const int output_argument = paired ? 4 : 3;
        SeedConfig config{50, 0, 1.0, 10000, 1000, 12, 5, 10, 0};
        for (int i = output_argument+1; i < argc;) {
            const std::string option = argv[i++];
            if (i >= argc) throw std::runtime_error("missing seed option value");
            const auto value = std::stoull(argv[i++]);
            if (option == "--seedSearchLmax") config.seedSearchLmax = value;
            else if (option == "--seedSearchStartLmax") config.seedSearchStartLmax = value;
            else if (option == "--seedMultimapNmax") config.seedMultimapNmax = value;
            else if (option == "--seedSplitMin") config.seedSplitMin = value;
            else if (option == "--seedMapMin") config.seedMapMin = value;
            else throw std::runtime_error("unsupported seed option: " + option);
        }
        std::ifstream first(argv[2]);
        std::ifstream second(paired ? argv[3] : "");
        std::ofstream output(argv[output_argument]);
        if (!first || (paired && !second) || !output)
            throw std::runtime_error("cannot open FASTQ or output");
        auto progress = StageProgress::fromEnvironment(
            std::string(argv[output_argument]) + ".timing.tsv",
            timedStageBit(TimedStage::Seed));
        auto reference = checked(loadGenomeReference(directory));
        SuffixArrayIndexBackend backend;
        SearchIndexRequest request{"nustar-sa", directory/"index.sa", {}};
        auto index = std::dynamic_pointer_cast<const SuffixArrayIndex>(
            checked(backend.load(request, reference)));
        if (!index) throw std::runtime_error("index type mismatch");
        StarSeedSearch search(index, reference, config);
        FastqRecord a, b;
        std::uint64_t ordinal = 0, piece_count = 0;
        double search_seconds = 0;
        while (readFastq(first, a)) {
            progress.beginRead();
            if (paired && !readFastq(second, b))
                throw std::runtime_error("FASTQ mates differ in length");
            const auto read = prepare(a, paired ? &b : nullptr, ordinal++);
            progress.beginStage(TimedStage::Seed);
            const auto started = std::chrono::steady_clock::now();
            const auto result = checked(search.search(read));
            search_seconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now()-started).count();
            progress.endStage(TimedStage::Seed);
            const auto& c = result.counters;
            output << "SEED_SUM\t" << a.name << '\t' << read.budget.legacy_Lread
                   << '\t' << result.splits.size() << '\t' << result.pieces.size()
                   << '\t' << c.all_alignments << '\t' << c.unique_alignments
                   << '\t' << c.multiple_alignments << '\t' << c.min_multimaps
                   << '\t' << c.min_multimap_length << '\t' << c.max_unique_length << '\n';
            for (std::size_t i = 0; i < result.splits.size(); ++i) {
                const auto& split = result.splits[i];
                output << "SEED_SPLIT\t" << a.name << '\t' << i << '\t'
                       << split.read.begin.value << '\t'
                       << split.read.end.value-split.read.begin.value << '\t'
                       << static_cast<unsigned>(split.mate) << '\n';
            }
            for (std::size_t i = 0; i < result.pieces.size(); ++i) {
                const auto& piece = result.pieces[i];
                output << "SEED_PC\t" << a.name << '\t' << i << '\t'
                       << piece.read_start << '\t' << piece.length << '\t'
                       << piece.direction << '\t' << piece.multiplicity << '\t'
                       << piece.sa_begin << '\t' << piece.sa_end << '\t'
                       << piece.fragment << '\n';
            }
            piece_count += result.pieces.size();
            progress.endRead();
        }
        if (paired && readFastq(second, b))
            throw std::runtime_error("FASTQ mates differ in length");
        output.close();
        progress.finish();
        std::cerr << "seed-trace reads=" << ordinal << " pieces=" << piece_count
                  << " search_seconds=" << search_seconds << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "seed-trace: " << error.what() << '\n';
        return 1;
    }
}
