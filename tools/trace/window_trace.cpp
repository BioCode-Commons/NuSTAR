#include "../../window/star_window.hpp"
#include "../../seed_search/star_search.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

using namespace nustar::star2711b;

namespace {

template<class T> T checked(Result<T> value) {
    if (std::holds_alternative<Failure>(value))
        throw std::runtime_error(std::get<Failure>(value).message);
    return std::get<T>(std::move(value));
}

void checked(Status value) {
    if (std::holds_alternative<Failure>(value))
        throw std::runtime_error(std::get<Failure>(value).message);
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

MappingRead prepare(const std::string& sequence, const std::string* mate2,
                    std::uint64_t ordinal) {
    MappingRead read{};
    read.fragment = FragmentId{0, ordinal};
    read.mates.push_back(PreparedMate{sequence.size(), sequence.size(), {},
        ReadRange{ReadOffset{0}, ReadOffset{sequence.size()}}, false});
    for (const auto base : sequence) read.encoded[0].push_back(encode(base));
    if (mate2) {
        read.encoded[0].push_back(BaseCode::MateSpacer);
        for (auto it = mate2->rbegin(); it != mate2->rend(); ++it) {
            const auto code = static_cast<std::uint8_t>(encode(*it));
            read.encoded[0].push_back(code < 4 ?
                static_cast<BaseCode>(code^3) : BaseCode::N);
        }
        const auto start = sequence.size()+1;
        read.mates.push_back(PreparedMate{mate2->size(), mate2->size(), {},
            ReadRange{ReadOffset{start}, ReadOffset{read.encoded[0].size()}}, true});
    }
    for (const auto base : read.encoded[0]) {
        const auto code = static_cast<std::uint8_t>(base);
        read.encoded[1].push_back(code < 4 ? static_cast<BaseCode>(code^3) : base);
    }
    read.encoded[2].assign(read.encoded[1].rbegin(), read.encoded[1].rend());
    read.budget = ReadBudget{10, read.encoded[0].size(),
        sequence.size()+(mate2 ? mate2->size() : 0),
        sequence.size()+(mate2 ? mate2->size() : 0)};
    return read;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 4 && argc != 5) {
        std::cerr << "Usage: window-trace NUSTAR_INDEX_DIR FIRST.fastq "
                     "[SECOND.fastq] OUTPUT.tsv\n";
        return 2;
    }
    try {
        auto reference = checked(loadGenomeReference(argv[1]));
        SuffixArrayIndexBackend backend;
        SearchIndexRequest request{"nustar-sa", std::filesystem::path(argv[1])/"index.sa", {}};
        auto index = checked(backend.load(request, reference));
        const SeedDeliveryContract delivery{SeedContract::Star2711bEquivalent,
            HitOrder::Star2711b, true, true, true, 1,
            std::string(star_seed_search_domain)};
        const SeedConfig seed_config{50, 0, 1.0, 10000, 1000, 12, 5, 10, 0};
        auto seed_plan = std::make_shared<SeedSearchPlan>(
            SeedSearchPlan{delivery, {}, seed_config, {}});
        StarSeedAlgorithmFactory seed_factory;
        auto finder = checked(seed_factory.bind(index, reference, seed_plan));
        auto builder = checked(makeStarWindowBuilder(reference));
        const auto& storage = reference->storage();
        const auto win_bits = std::min<std::uint32_t>(16, storage.chromosome_bin_bits);
        const WindowConfig window_config{50, win_bits,
            storage.chromosome_bin_bits-win_bits,
            (storage.forward.size()+reference->synthetic().size())/(1ULL<<win_bits)+1,
            9, 4, 50, 10, 10000, 0, 0};
        const CompatibilityLimits limits{650, 1024, 20, 10000, 255};
        const WorkspaceBudget budget{std::numeric_limits<std::uint64_t>::max(), 128};
        auto seed_workspace = checked(finder->createWorkspace(limits, budget));
        auto window_workspace = checked(builder->createWorkspace(
            window_config, reference->view(), limits, budget));
        const bool paired = argc == 5;
        std::ifstream input(argv[2]);
        std::ifstream second(paired ? argv[3] : "");
        std::ofstream output(argv[paired ? 4 : 3]);
        if (!input || (paired && !second) || !output)
            throw std::runtime_error("cannot open FASTQ or output");
        std::string name, sequence, plus, quality;
        std::string mate_name, mate_sequence, mate_plus, mate_quality;
        std::uint64_t ordinal = 0;
        while (std::getline(input, name)) {
            if (!std::getline(input, sequence) || !std::getline(input, plus) ||
                !std::getline(input, quality) || name.empty() || name[0] != '@' ||
                plus.empty() || plus[0] != '+' || quality.size() != sequence.size())
                throw std::runtime_error("invalid FASTQ record");
            if (paired && (!std::getline(second, mate_name) ||
                !std::getline(second, mate_sequence) ||
                !std::getline(second, mate_plus) ||
                !std::getline(second, mate_quality) || mate_name.empty() ||
                mate_name[0] != '@' || mate_plus.empty() ||
                mate_plus[0] != '+' ||
                mate_quality.size() != mate_sequence.size()))
                throw std::runtime_error("invalid mate 2 FASTQ record");
            const auto name_end = name.find_first_of(" \t");
            if (name_end != std::string::npos) name.resize(name_end);
            auto read = prepare(sequence, paired ? &mate_sequence : nullptr,
                                ordinal++);
            SeedBatch seeds{};
            checked(finder->find(read, *seed_workspace, seeds));
            WindowBatch windows{};
            checked(builder->build(read, seeds, reference->view(), window_config,
                                   *window_workspace, windows));
            output << "STITCH_READ\t" << name << '\t' << read.budget.legacy_Lread
                   << "\t10\n";
            for (const auto& window : windows.windows) {
                if (!window.alive || window.seeds.empty()) continue;
                output << "STITCH_WIN\t" << name << '\t' << window.id.value
                       << '\t' << window.chromosome.value << '\t'
                       << static_cast<unsigned>(window.orientation) << '\t'
                       << static_cast<unsigned>(window.orientation) << '\t'
                       << read.budget.legacy_Lread << "\t10\t" << window.seeds.size()
                       << '\t' << (window.last_anchor ? *window.last_anchor :
                           std::numeric_limits<std::uint64_t>::max()) << '\n';
                for (std::size_t i = 0; i < window.seeds.size(); ++i) {
                    const auto& seed = window.seeds[i];
                    output << "STITCH_WA\t" << name << '\t' << window.id.value
                           << '\t' << i << '\t' << seed.read_start.value << '\t'
                           << seed.genome_start.value << '\t' << seed.length << '\t'
                           << seed.repeat_count << '\t' << seed.anchor << '\t'
                           << static_cast<unsigned>(seed.mate) << '\t'
                           << (seed.sjdb ? seed.sjdb->value :
                               std::numeric_limits<std::uint64_t>::max()) << '\n';
                }
            }
        }
        if (paired && std::getline(second, mate_name))
            throw std::runtime_error("mate 2 contains extra records");
        std::cerr << "window-trace reads=" << ordinal << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "window-trace: " << error.what() << '\n';
        return 1;
    }
}
