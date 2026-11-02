#include "gtf_loader.hpp"
#include <algorithm>
#include <charconv>
#include <fstream>
#include <map>
#include <tuple>
#include <unordered_map>

namespace nustar::star2711b {

namespace {
Failure gtfError(ErrorKind kind, std::string message) {
    return Failure{kind, Stage::Gtf, std::nullopt, 1, std::move(message)};
}

std::string_view trim(std::string_view value) {
    const auto first = value.find_first_not_of(" \t\r\n\"");
    if (first == std::string_view::npos) return {};
    const auto last = value.find_last_not_of(" \t\r\n\"");
    return value.substr(first, last - first + 1);
}

std::map<std::string, std::string> parseAttributes(std::string_view text) {
    std::map<std::string, std::string> attrs;
    while (!text.empty()) {
        const auto semi = text.find(';');
        auto field = trim(text.substr(0, semi));
        const auto split = field.find_first_of(" \t=");
        if (split != std::string_view::npos) {
            const auto key = trim(field.substr(0, split));
            const auto value = trim(field.substr(split + 1));
            if (!key.empty()) attrs.emplace(key, value);
        }
        if (semi == std::string_view::npos) break;
        text.remove_prefix(semi + 1);
    }
    return attrs;
}

std::string choose(const std::map<std::string, std::string>& attrs,
                   const std::vector<std::string>& names) {
    for (const auto& name : names) {
        if (const auto it = attrs.find(name); it != attrs.end() && !it->second.empty())
            return it->second;
    }
    return {};
}

std::uint64_t number(std::string_view text) {
    std::uint64_t value = 0;
    const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
    return error == std::errc{} && end == text.data() + text.size() ? value : 0;
}

} // namespace

Result<AnnotationCatalog> GtfAnnotationLoader::load(
    const AnnotationConfig& config, const ReferenceView& reference) {
    AnnotationCatalog catalog{};
    if (!config.sjdbGTFfile) return catalog;
    if (config.sjdbOverhang == 0)
        return gtfError(ErrorKind::Config, "GTF requires positive sjdbOverhang");
    std::ifstream input(*config.sjdbGTFfile);
    if (!input) return gtfError(ErrorKind::Io, "cannot open GTF: " + config.sjdbGTFfile->string());
    catalog.active = true;
    std::unordered_map<std::string, std::uint64_t> chromosome_ids, gene_ids, transcript_ids;
    for (const auto& chr : reference.metadata->chromosomes) chromosome_ids.emplace(chr.name, chr.id.value);
    std::string line;
    std::uint64_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (line.empty() || line[0] == '#') continue;
        std::array<std::string_view, 9> columns{};
        std::string_view row(line);
        for (std::size_t i = 0; i < 8; ++i) {
            const auto tab = row.find('\t');
            if (tab == std::string_view::npos) break;
            columns[i] = row.substr(0, tab);
            row.remove_prefix(tab + 1);
            if (i == 7) columns[8] = row;
        }
        if (columns[2] != config.sjdbGTFfeatureExon) continue;
        if (columns[8].empty()) {
            catalog.diagnostics.push_back(Diagnostic{Diagnostic::Level::Warning, Stage::Gtf,
                *config.sjdbGTFfile, line_number, "gtf_columns", "invalid exon row"});
            continue;
        }
        const auto chr_name = config.sjdbGTFchrPrefix + std::string(columns[0]);
        const auto chr_it = chromosome_ids.find(chr_name);
        if (chr_it == chromosome_ids.end()) {
            catalog.diagnostics.push_back(Diagnostic{Diagnostic::Level::Warning, Stage::Gtf,
                *config.sjdbGTFfile, line_number, "unknown_chromosome", chr_name});
            continue;
        }
        const auto first = number(columns[3]), last = number(columns[4]);
        const auto& chr = reference.metadata->chromosomes[chr_it->second];
        if (first == 0 || last < first || last > chr.length) {
            catalog.diagnostics.push_back(Diagnostic{Diagnostic::Level::Warning, Stage::Gtf,
                *config.sjdbGTFfile, line_number, "exon_coordinates", "exon outside chromosome"});
            continue;
        }
        const auto attrs = parseAttributes(columns[8]);
        auto transcript = choose(attrs, {config.sjdbGTFtagExonParentTranscript});
        if (transcript.empty()) {
            transcript = "tr_" + chr_name + "_" + std::to_string(first) + "_" +
                         std::to_string(last) + "_" + std::to_string(catalog.exons.size());
            catalog.diagnostics.push_back(Diagnostic{Diagnostic::Level::Warning, Stage::Gtf,
                *config.sjdbGTFfile, line_number, "missing_transcript_id", transcript});
        }
        auto gene = choose(attrs, {config.sjdbGTFtagExonParentGene});
        if (gene.empty()) {
            gene = "MissingGeneID";
            catalog.diagnostics.push_back(Diagnostic{Diagnostic::Level::Warning, Stage::Gtf,
                *config.sjdbGTFfile, line_number, "missing_gene_id", gene});
        }
        auto gene_name = choose(attrs, config.sjdbGTFtagExonParentGeneName);
        if (gene_name.empty()) gene_name = gene;
        auto gene_type = choose(attrs, config.sjdbGTFtagExonParentGeneType);
        if (gene_type.empty()) gene_type = "MissingGeneType";
        std::uint64_t gene_id = 0, transcript_id = 0;
        if (const auto [it, inserted] = gene_ids.emplace(gene, catalog.genes.size()); inserted) {
            gene_id = it->second;
            catalog.genes.push_back(GeneModel{GeneId{gene_id}, gene, gene_name, gene_type});
        } else gene_id = it->second;
        if (const auto [it, inserted] = transcript_ids.emplace(transcript, catalog.transcripts.size()); inserted) {
            transcript_id = it->second;
            const auto strand = columns[6] == "+" ? Strand::Plus :
                                columns[6] == "-" ? Strand::Minus : Strand::Unknown;
            catalog.transcripts.push_back(TranscriptModel{TranscriptId{transcript_id},
                transcript, strand, {}});
        } else transcript_id = it->second;
        const auto index = catalog.exons.size();
        catalog.exons.push_back(ExonModel{TranscriptId{transcript_id}, GeneId{gene_id},
            chr.id, {GenomeOffset{chr.start.value + first - 1},
                     GenomeOffset{chr.start.value + last}}, line_number, index});
        catalog.transcripts[transcript_id].exon_indices.push_back(index);
    }
    if (!input.eof()) return gtfError(ErrorKind::Io, "error reading GTF");
    return catalog;
}

Result<std::vector<JunctionProposal>> GtfAnnotationLoader::deriveJunctions(
    const AnnotationCatalog& catalog, const ReferenceView&) {
    std::vector<JunctionProposal> proposals;
    std::map<std::tuple<std::uint64_t,std::uint64_t,std::uint64_t,std::uint8_t>, std::size_t> known;
    for (const auto& transcript : catalog.transcripts) {
        auto indices = transcript.exon_indices;
        std::stable_sort(indices.begin(), indices.end(), [&](auto a, auto b) {
            const auto& x = catalog.exons[a];
            const auto& y = catalog.exons[b];
            return std::tie(x.interval.begin.value, x.interval.end.value) <
                   std::tie(y.interval.begin.value, y.interval.end.value);
        });
        for (std::size_t i = 1; i < indices.size(); ++i) {
            const auto& previous = catalog.exons[indices[i - 1]];
            const auto& current = catalog.exons[indices[i]];
            if (previous.chromosome.value != current.chromosome.value ||
                current.interval.begin.value <= previous.interval.end.value) continue;
            const auto key = std::make_tuple(current.chromosome.value,
                previous.interval.end.value, current.interval.begin.value,
                static_cast<std::uint8_t>(transcript.strand));
            const auto [it, inserted] = known.emplace(key, proposals.size());
            if (inserted) proposals.push_back(JunctionProposal{current.chromosome,
                {previous.interval.end, current.interval.begin}, transcript.strand,
                JunctionOrigin::Gtf, 20, current.encounter_order, {current.gene}});
            else proposals[it->second].genes.push_back(current.gene);
        }
    }
    for (auto& proposal : proposals) {
        std::sort(proposal.genes.begin(), proposal.genes.end(),
                  [](GeneId a, GeneId b) { return a.value < b.value; });
        proposal.genes.erase(std::unique(proposal.genes.begin(), proposal.genes.end(),
            [](GeneId a, GeneId b) { return a.value == b.value; }), proposal.genes.end());
    }
    std::stable_sort(proposals.begin(), proposals.end(), [](const auto& a, const auto& b) {
        return std::tie(a.intron.begin.value, a.intron.end.value, a.strand, a.source_order) <
               std::tie(b.intron.begin.value, b.intron.end.value, b.strand, b.source_order);
    });
    return proposals;
}

} // namespace nustar::star2711b
