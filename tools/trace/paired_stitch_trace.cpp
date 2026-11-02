#include "../../stitching/paired_stitch.hpp"
#include "../../reference/genome.hpp"
#include "../../profiling/stage_progress.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>

using namespace nustar::star2711b;

namespace {

template<class T> T checked(Result<T> value) {
    if (std::holds_alternative<Failure>(value))
        throw std::runtime_error(std::get<Failure>(value).message);
    return std::get<T>(std::move(value));
}

std::vector<std::string> fields(const std::string& line) {
    std::vector<std::string> result;
    std::size_t begin = 0;
    while (begin <= line.size()) {
        const auto end = line.find('\t', begin);
        result.push_back(line.substr(begin, end == std::string::npos ? end : end-begin));
        if (end == std::string::npos) break;
        begin = end+1;
    }
    return result;
}

std::uint64_t number(const std::vector<std::string>& values, std::size_t index) {
    if (index >= values.size()) throw std::runtime_error("short STAR stitch trace row");
    return std::stoull(values[index]);
}

BaseCode base(char letter) {
    switch (letter) {
    case 'A': case 'a': return BaseCode::A;
    case 'C': case 'c': return BaseCode::C;
    case 'G': case 'g': return BaseCode::G;
    case 'T': case 't': return BaseCode::T;
    default: return BaseCode::N;
    }
}

std::unordered_map<std::string, std::string> loadFastq(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open FASTQ");
    std::unordered_map<std::string, std::string> reads;
    std::string name, sequence, plus, quality;
    while (std::getline(input, name)) {
        if (!std::getline(input, sequence) || !std::getline(input, plus) ||
            !std::getline(input, quality) || name.empty() || name[0] != '@' ||
            plus.empty() || plus[0] != '+' || quality.size() != sequence.size())
            throw std::runtime_error("malformed FASTQ");
        const auto name_end = name.find_first_of(" \t");
        if (name_end != std::string::npos) name.resize(name_end);
        reads.emplace(std::move(name), std::move(sequence));
    }
    return reads;
}

PairedMappingRead makeRead(const std::string& first,
                           const std::string& second,
                           std::uint64_t mismatch_limit, bool revert) {
    MappingRead read{};
    read.fragment = FragmentId{0, 0};
    read.mates.push_back(PreparedMate{first.size(), first.size(), {},
        ReadRange{ReadOffset{0}, ReadOffset{first.size()}}, false});
    const auto total = first.size()+second.size()+1;
    read.budget = ReadBudget{mismatch_limit, total, first.size()+second.size(),
                             first.size()+second.size()};
    read.revert_strand = revert;
    for (const auto letter : first) read.encoded[0].push_back(base(letter));
    read.encoded[0].push_back(BaseCode::MateSpacer);
    for (auto it = second.rbegin(); it != second.rend(); ++it) {
        const auto value = static_cast<std::uint8_t>(base(*it));
        read.encoded[0].push_back(value < 4 ?
            static_cast<BaseCode>(value^3) : BaseCode::N);
    }
    read.mates.push_back(PreparedMate{second.size(), second.size(), {},
        ReadRange{ReadOffset{first.size()+1}, ReadOffset{total}}, true});
    for (const auto code : read.encoded[0]) {
        const auto value = static_cast<std::uint8_t>(code);
        read.encoded[1].push_back(value < 4 ? static_cast<BaseCode>(value^3) : code);
    }
    read.encoded[2].assign(read.encoded[1].rbegin(), read.encoded[1].rend());
    return PairedMappingRead::fromPrepared(std::move(read));
}

PairedStitchConfig starConfig() {
    PairedStitchConfig result{};
    auto& c = result.common;
    c.alignIntronMin = 21;
    c.alignIntronMax = 0;
    c.alignSJoverhangMin = 5;
    c.alignSJDBoverhangMin = 3;
    c.alignSJstitchMismatchNmax = {0, -1, 0, 0};
    c.alignSplicedMateMapLmin = 0;
    c.alignSplicedMateMapLminOverLmate = 0.66;
    c.alignTranscriptsPerWindowNmax = 100;
    c.alignTranscriptsPerReadNmax = 10000;
    c.alignEndsType = "Local";
    c.alignSoftClipAtReferenceEnds = true;
    c.insertion_flush_right = false;
    c.outFilterIntronMotifs = "None";
    c.outFilterIntronStrands = "RemoveInconsistentStrands";
    c.outSAMstrandField = "None";
    c.outFilterMismatchNoverLmax = 0.3;
    c.outFilterMultimapScoreRange = 1;
    result.force_end_extension = {{{false, false}, {false, false}}};
    result.alignMatesGapMax = 0;
    result.protrusion_bases_max = 0;
    result.protrusion_concordant_pair = false;
    return result;
}

ScoreConfig starScore() {
    ScoreConfig result{};
    result.scoreGap = 0;
    result.scoreGapNoncan = -8;
    result.scoreGapGCAG = -4;
    result.scoreGapATAC = -8;
    result.scoreDelOpen = -2;
    result.scoreDelBase = -2;
    result.scoreInsOpen = -2;
    result.scoreInsBase = -2;
    result.scoreStitchSJshift = 1;
    result.sjdbScore = 2;
    result.scoreGenomicLengthLog2scale = -0.25;
    return result;
}

int boundaryCode(const AlignmentBoundary& boundary) {
    switch (boundary.kind) {
    case BoundaryKind::Insertion: return -2;
    case BoundaryKind::Deletion: return -1;
    case BoundaryKind::MateGap: return -3;
    case BoundaryKind::Intron: return static_cast<int>(boundary.motif);
    }
    throw std::runtime_error("unknown boundary kind");
}

void emit(std::ostream& output, const std::string& name, std::uint64_t window,
          const CandidateSet& result, std::int64_t nanoseconds) {
    output << "STITCH_OUT\t" << name << '\t' << window << '\t'
           << result.candidates.size() << '\t' << nanoseconds << '\n';
    for (std::size_t i = 0; i < result.candidates.size(); ++i) {
        const auto& c = result.candidates[i];
        const auto& m = c.metrics;
        output << "STITCH_TR\t" << name << '\t' << window << '\t' << i << '\t'
               << m.score << '\t' << c.read_start.value << '\t'
               << c.original_read_start.value << '\t' << c.read_length << '\t'
               << c.genome_start.value << '\t' << c.genome_length << '\t'
               << c.blocks.size() << '\t' << m.matches << '\t' << m.mismatches
               << '\t' << m.mapped_length << '\t' << m.gap_count << '\t'
               << m.gap_bases << '\t' << m.deletion_count << '\t'
               << m.deletion_bases << '\t' << m.insertion_count << '\t'
               << m.insertion_bases << '\t' << m.anchor_seeds << '\n';
        output << "STITCH_META\t" << name << '\t' << window << '\t' << i
               << '\t' << m.unique_seeds << '\t' << m.extension_length
               << '\t' << m.intron_motif_counts[0] << '\t' << m.intron_motif_counts[1]
               << '\t' << m.intron_motif_counts[2] << '\t'
               << unsigned(m.motif_strand) << '\t' << unsigned(m.has_splice)
               << '\t' << (c.single_mate ? static_cast<int>(*c.single_mate) : -1)
               << '\n';
        for (std::size_t exon = 0; exon < c.blocks.size(); ++exon) {
            const auto& b = c.blocks[exon];
            output << "STITCH_EXON\t" << name << '\t' << window << '\t' << i
                   << '\t' << exon << '\t' << b.read_start.value << '\t'
                   << b.genome_start.value << '\t' << b.length << '\t'
                   << static_cast<unsigned>(b.mate) << '\t'
                   << (b.sjdb ? b.sjdb->value : std::numeric_limits<std::uint64_t>::max())
                   << '\n';
            if (exon < c.boundaries.size()) {
                const auto& j = c.boundaries[exon];
                output << "STITCH_JUNCTION\t" << name << '\t' << window << '\t'
                       << i << '\t' << exon << '\t' << boundaryCode(j) << '\t'
                       << unsigned(j.annotation) << '\t'
                       << unsigned(j.strand) << '\t'
                       << j.shift_left << '\t' << j.shift_right << '\n';
            }
        }
    }
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "Usage: paired-stitch-trace INDEX_DIR MATE1.fastq "
                     "MATE2.fastq STAR_TRACE.stderr OUTPUT.tsv\n";
        return 2;
    }
    try {
        auto progress = StageProgress::fromEnvironment(
            std::string(argv[5]) + ".timing.tsv",
            timedStageBit(TimedStage::Stitch));
        const auto reads = loadFastq(argv[2]);
        const auto mates = loadFastq(argv[3]);
        const auto reference = checked(loadGenomeReference(argv[1]));
        const auto view = reference->view();
        const auto stitcher = checked(makeStarPairedStitcher(reference));
        const auto config = starConfig();
        const auto score = starScore();
        const CompatibilityLimits limits{650, 1024, 20, 10000, 255};
        const WorkspaceBudget budget{std::numeric_limits<std::uint64_t>::max(), 4096};
        auto workspace = checked(stitcher->createWorkspace(config, score, view, limits, budget));
        std::ifstream trace(argv[4]);
        std::ofstream output(argv[5]);
        if (!trace || !output) throw std::runtime_error("cannot open stitch trace files");
        std::string line, current_name;
        AlignmentWindow window{};
        std::vector<AlignmentWindow> current_windows;
        std::uint64_t mismatch_limit = 0, read_length = 0, strand = 0, original_strand = 0;
        std::size_t expected_seeds = 0;
        std::uint64_t processed = 0, total_candidates = 0;
        double seconds = 0;
        while (std::getline(trace, line)) {
            if (line.rfind("STITCH_READ\t", 0) == 0) {
                progress.beginRead();
                const auto f = fields(line);
                if (f.size() != 4 || !current_windows.empty())
                    throw std::runtime_error("invalid STITCH_READ row");
                current_name = f[1];
                read_length = number(f, 2);
                mismatch_limit = number(f, 3);
                strand = original_strand = 0;
            } else if (line.rfind("STITCH_WIN\t", 0) == 0) {
                const auto f = fields(line);
                if (f.size() != 10 || f[1] != current_name ||
                    number(f, 6) != read_length || number(f, 7) != mismatch_limit)
                    throw std::runtime_error("invalid STITCH_WIN row");
                window = AlignmentWindow{};
                window.id = WindowId{number(f, 2)};
                window.chromosome = ChromosomeId{number(f, 3)};
                strand = number(f, 4);
                original_strand = number(f, 5);
                window.orientation = static_cast<Orientation>(strand);
                read_length = number(f, 6);
                mismatch_limit = number(f, 7);
                expected_seeds = number(f, 8);
                const auto last_anchor = number(f, 9);
                if (last_anchor != std::numeric_limits<std::uint64_t>::max())
                    window.last_anchor = last_anchor;
                window.alive = true;
            } else if (line.rfind("STITCH_WA\t", 0) == 0) {
                const auto f = fields(line);
                if (f.size() != 11 || f[1] != current_name ||
                    number(f, 2) != window.id.value ||
                    number(f, 3) != window.seeds.size())
                    throw std::runtime_error("invalid STITCH_WA row");
                const auto sjdb = number(f, 10);
                const auto anchor = number(f, 8);
                window.seeds.push_back(WindowSeed{ReadOffset{number(f, 4)},
                    GenomeOffset{number(f, 5)}, number(f, 6), number(f, 7),
                    anchor != 0, static_cast<MateId>(number(f, 9)),
                    sjdb == std::numeric_limits<std::uint64_t>::max() ?
                        std::optional<JunctionId>{} :
                        std::optional<JunctionId>{JunctionId{sjdb}}});
            } else if (line.rfind("STITCH_OUT\t", 0) == 0) {
                if (window.seeds.size() != expected_seeds)
                    throw std::runtime_error("STAR stitch window seed count differs");
                const auto found = reads.find(current_name);
                const auto mate = mates.find(current_name);
                if (found == reads.end() || mate == mates.end() ||
                    found->second.size()+mate->second.size()+1 != read_length)
                    throw std::runtime_error("STAR paired read not present in FASTQ");
                const auto read = makeRead(found->second, mate->second, mismatch_limit,
                                           original_strand != strand);
                WindowBatch batch{};
                batch.generation = view.generation;
                batch.windows = current_windows;
                batch.windows.push_back(window);
                CandidateSet complete;
                progress.beginStage(TimedStage::Stitch);
                const auto started = std::chrono::steady_clock::now();
                const auto status = stitcher->stitchPair(read, batch, view, config,
                    score, *workspace, complete);
                const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now()-started).count();
                progress.endStage(TimedStage::Stitch);
                if (std::holds_alternative<Failure>(status))
                    throw std::runtime_error(std::get<Failure>(status).message);
                CandidateSet result{};
                for (const auto& wc : complete.windows)
                    if (wc.window.value == window.id.value)
                        for (const auto id : wc.ordered)
                            result.candidates.push_back(complete.candidates[id.value]);
                emit(output, current_name, window.id.value, result, elapsed);
                current_windows.push_back(window);
                ++processed;
                total_candidates += result.candidates.size();
                seconds += double(elapsed)/1e9;
            } else if (line.rfind("STITCH_BEST\t", 0) == 0) {
                const auto f = fields(line);
                if (f.size() != 8 || f[1] != current_name)
                    throw std::runtime_error("invalid STITCH_BEST row");
                const auto found = reads.find(current_name);
                const auto mate = mates.find(current_name);
                if (found == reads.end() || mate == mates.end() ||
                    found->second.size()+mate->second.size()+1 != read_length)
                    throw std::runtime_error("STAR paired read not present in FASTQ");
                const auto read = makeRead(found->second, mate->second, mismatch_limit,
                                           original_strand != strand);
                WindowBatch batch{};
                batch.generation = view.generation;
                batch.windows = std::move(current_windows);
                CandidateSet result;
                const auto status = stitcher->stitchPair(read, batch, view, config,
                    score, *workspace, result);
                if (std::holds_alternative<Failure>(status))
                    throw std::runtime_error(std::get<Failure>(status).message);
                if (result.best) {
                    const auto& best = result.candidates[result.best->value];
                    output << "STITCH_BEST\t" << current_name << '\t'
                           << best.metrics.score << '\t' << best.genome_start.value
                           << '\t' << best.genome_length << '\t'
                           << best.chromosome.value << '\t'
                           << static_cast<unsigned>(best.orientation) << '\t'
                           << result.windows.size() << '\n';
                } else {
                    output << "STITCH_BEST\t" << current_name
                           << "\t0\t0\t0\t0\t0\t0\n";
                }
                current_windows.clear();
                progress.endRead();
            }
        }
        progress.finish();
        std::cerr << "paired-stitch-trace windows=" << processed
                  << " candidates=" << total_candidates
                  << " stitch_seconds=" << seconds << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "paired-stitch-trace: " << error.what() << '\n';
        return 1;
    }
}
