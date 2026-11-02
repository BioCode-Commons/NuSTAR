#include "table_writer.hpp"
#include <algorithm>
#include <fstream>
#include <tuple>

namespace nustar::star2711b {

namespace {
Failure gtfError(ErrorKind kind, std::string message) {
    return Failure{kind, Stage::Gtf, std::nullopt, 1, std::move(message)};
}

char strandChar(Strand strand) {
    return strand == Strand::Plus ? '+' : strand == Strand::Minus ? '-' : '.';
}

} // namespace

Status GtfTableWriter::writeLegacyTables(const AnnotationCatalog& catalog,
    std::span<const JunctionProposal> proposals, const ReferenceMetadata& metadata,
    const std::filesystem::path& directory) {
    std::error_code error;
    std::filesystem::create_directories(directory, error);
    if (error) return gtfError(ErrorKind::Io, "cannot create annotation directory");
    std::ofstream genes(directory / "geneInfo.tab");
    std::ofstream exons_by_gene(directory / "exonGeTrInfo.tab");
    std::ofstream transcripts(directory / "transcriptInfo.tab");
    std::ofstream exons(directory / "exonInfo.tab");
    std::ofstream junctions(directory / "sjdbList.fromGTF.out.tab");
    if (!genes || !exons_by_gene || !transcripts || !exons || !junctions)
        return gtfError(ErrorKind::Io, "cannot create GTF auxiliary tables");
    genes << catalog.genes.size() << '\n';
    for (const auto& gene : catalog.genes)
        genes << gene.external_id << '\t' << gene.name << '\t' << gene.biotype << '\n';
    std::vector<const ExonModel*> sorted_exons;
    for (const auto& exon : catalog.exons) sorted_exons.push_back(&exon);
    std::stable_sort(sorted_exons.begin(), sorted_exons.end(), [&](auto a, auto b) {
        // STAR funCompareArrays<uint64,5> compares every emitted column.
        return std::tuple(a->interval.begin.value, a->interval.end.value,
                         catalog.transcripts[a->transcript.value].strand, a->gene.value, a->transcript.value) <
               std::tuple(b->interval.begin.value, b->interval.end.value,
                         catalog.transcripts[b->transcript.value].strand, b->gene.value, b->transcript.value);
    });
    exons_by_gene << catalog.exons.size() << '\n';
    for (const auto* exon : sorted_exons)
        exons_by_gene << exon->interval.begin.value << '\t' << exon->interval.end.value - 1 << '\t'
            << static_cast<int>(catalog.transcripts[exon->transcript.value].strand) << '\t'
            << exon->gene.value << '\t' << exon->transcript.value << '\n';
    transcripts << catalog.transcripts.size() << '\n';
    exons << catalog.exons.size() << '\n';
    std::vector<const TranscriptModel*> sorted_transcripts;
    sorted_transcripts.reserve(catalog.transcripts.size());
    for (const auto& transcript : catalog.transcripts)
        if (!transcript.exon_indices.empty()) sorted_transcripts.push_back(&transcript);
    auto bounds = [&](const TranscriptModel& transcript) {
        std::uint64_t first = UINT64_MAX, last = 0;
        for (auto index : transcript.exon_indices) {
            first = std::min(first, catalog.exons[index].interval.begin.value);
            last = std::max(last, catalog.exons[index].interval.end.value);
        }
        return std::pair{first, last};
    };
    std::stable_sort(sorted_transcripts.begin(), sorted_transcripts.end(),
        [&](auto a, auto b) {
            const auto [as, ae] = bounds(*a);
            const auto [bs, be] = bounds(*b);
            return std::tie(as, ae, a->id.value) < std::tie(bs, be, b->id.value);
        });
    std::uint64_t exon_ordinal = 0;
    std::uint64_t previous_max_end = 0;
    for (const auto* transcript_ptr : sorted_transcripts) {
        const auto& transcript = *transcript_ptr;
        auto indices = transcript.exon_indices;
        std::stable_sort(indices.begin(), indices.end(), [&](auto a, auto b) {
            return catalog.exons[a].interval.begin.value < catalog.exons[b].interval.begin.value;
        });
        const auto first = catalog.exons[indices.front()].interval.begin.value;
        std::uint64_t last = first, spliced = 0;
        for (auto index : indices) last = std::max(last, catalog.exons[index].interval.end.value);
        if (exon_ordinal == 0) previous_max_end = last - 1;
        transcripts << transcript.external_id << '\t' << first << '\t' << last - 1 << '\t'
            << previous_max_end << '\t' << static_cast<int>(transcript.strand) << '\t'
            << indices.size() << '\t' << exon_ordinal << '\t'
            << catalog.exons[indices.back()].gene.value << '\n';
        for (auto index : indices) {
            const auto& exon = catalog.exons[index];
            exons << exon.interval.begin.value - first << '\t'
                << exon.interval.end.value - 1 - first << '\t' << spliced << '\n';
            spliced += exon.interval.end.value - exon.interval.begin.value;
            ++exon_ordinal;
        }
        previous_max_end = std::max(previous_max_end, last - 1);
    }
    for (const auto& proposal : proposals) {
        const auto& chr = metadata.chromosomes.at(proposal.chromosome.value);
        junctions << chr.name << '\t' << proposal.intron.begin.value - chr.start.value + 1
                  << '\t' << proposal.intron.end.value - chr.start.value << '\t'
                  << strandChar(proposal.strand) << '\t';
        for (std::size_t i = 0; i < proposal.genes.size(); ++i) {
            if (i) junctions << ',';
            junctions << proposal.genes[i].value + 1;
        }
        junctions << '\n';
    }
    if (!genes || !exons_by_gene || !transcripts || !exons || !junctions)
        return gtfError(ErrorKind::Io, "error writing GTF auxiliary tables");
    return std::monostate{};
}

namespace {
const Chromosome* chromosomeAt(const ReferenceMetadata& metadata, GenomeOffset position) {
    for (const auto& chr : metadata.chromosomes) {
        if (position.value >= chr.start.value && position.value < chr.start.value + chr.length)
            return &chr;
    }
    return nullptr;
}

} // namespace

Status writeJunctionDatabaseTables(const ReferenceView& reference,
    const std::filesystem::path& directory) {
    const auto info_path = directory / "sjdbInfo.txt";
    const auto list_path = directory / "sjdbList.out.tab";
    std::ofstream info(info_path), list(list_path);
    if (!info || !list)
        return Failure{ErrorKind::Io, Stage::Sjdb, {}, 1, "cannot write junction tables"};
    info << reference.junctions.entries.size() << '\t' << reference.junctions.overhang << '\n';
    for (const auto& j : reference.junctions.entries) {
        info << j.lookup_intron.begin.value << '\t' << j.lookup_intron.end.value - 1 << '\t'
             << static_cast<int>(j.motif) << '\t' << static_cast<int>(j.shift_left) << '\t'
             << static_cast<int>(j.shift_right) << '\t' << static_cast<int>(j.strand) << '\n';
        const auto* chr = chromosomeAt(*reference.metadata, j.reported_intron.begin);
        if (!chr) return Failure{ErrorKind::Internal, Stage::Sjdb, {}, 1,
                                      "junction lost its chromosome"};
        list << chr->name << '\t'
             << j.reported_intron.begin.value - chr->start.value + 1 << '\t'
             << j.reported_intron.end.value - chr->start.value << '\t'
             << (j.strand == Strand::Plus ? '+' : j.strand == Strand::Minus ? '-' : '.') << '\n';
    }
    info.close(); list.close();
    if (!info || !list) return Failure{ErrorKind::Io, Stage::Sjdb, {}, 1,
                                             "error writing junction tables"};
    return std::monostate{};
}

} // namespace nustar::star2711b
