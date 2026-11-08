#include "interface.hpp"

#include <ctime>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>

namespace nustar::star2711b {

Result<StatsDelta> deriveStatistics(const MappingRead& read,
    const CandidateSet& candidates, const SelectionResult& selection,
    const FilterDecision& decision, const OutputDisposition& disposition) {
    StatsDelta delta{};
    delta.reads = 1;
    for (const auto& mate : read.mates)
        delta.read_bases += mate.clipped_length; // STAR counts readLength after clipping.

    switch (decision.reason) {
    case UnmappedReason::None: break;
    case UnmappedReason::NoGoodWindow: ++delta.unmapped_other; break;
    case UnmappedReason::TooShort: ++delta.unmapped_short; break;
    case UnmappedReason::TooManyMismatches: ++delta.unmapped_mismatch; break;
    case UnmappedReason::TooManyLoci: ++delta.unmapped_multi; break;
    }
    if (disposition.count_as_unmapped_all) ++delta.unmapped_all;

    if (decision.mapping == MappingClass::Multiple) {
        ++delta.mapped_multiple;
    } else if (decision.mapping == MappingClass::Unique) {
        ++delta.mapped_unique;
        // Unique mapping implies a selected candidate and a nonempty mapping read.
        const auto& alignment = candidates.candidates[
            selection.ordered.front().candidate.value];
        const auto& metrics = alignment.metrics;
        delta.mismatches = metrics.mismatches;
        delta.insertion_count = metrics.insertion_count;
        delta.deletion_count = metrics.deletion_count;
        delta.insertion_bases = metrics.insertion_bases;
        delta.deletion_bases = metrics.deletion_bases;
        for (const auto& block : alignment.blocks)
            delta.mapped_bases += block.length;
        delta.mapped_portion = double(delta.mapped_bases) /
                               double(read.budget.legacy_Lread);
        for (const auto& boundary : alignment.boundaries) {
            if (boundary.kind != BoundaryKind::Intron) continue;
            ++delta.splices_by_motif[static_cast<std::size_t>(boundary.motif)];
            if (boundary.annotation == 1) ++delta.annotated_splices;
        }
    }
    return delta;
}

void addStatistics(StatsDelta& total, const StatsDelta& delta) noexcept {
    total.reads += delta.reads;
    total.read_bases += delta.read_bases;
    total.mapped_unique += delta.mapped_unique;
    total.mapped_multiple += delta.mapped_multiple;
    total.mapped_bases += delta.mapped_bases;
    total.mismatches += delta.mismatches;
    total.insertion_count += delta.insertion_count;
    total.deletion_count += delta.deletion_count;
    total.insertion_bases += delta.insertion_bases;
    total.deletion_bases += delta.deletion_bases;
    total.mapped_portion += delta.mapped_portion;
    for (std::size_t i = 0; i < total.splices_by_motif.size(); ++i)
        total.splices_by_motif[i] += delta.splices_by_motif[i];
    total.annotated_splices += delta.annotated_splices;
    total.unmapped_other += delta.unmapped_other;
    total.unmapped_short += delta.unmapped_short;
    total.unmapped_mismatch += delta.unmapped_mismatch;
    total.unmapped_multi += delta.unmapped_multi;
    total.unmapped_all += delta.unmapped_all;
}

QualityHistogramDelta deriveQualityHistogram(
    std::span<const std::string_view> quality_by_mate) {
    QualityHistogramDelta result{};
    for (std::size_t mate = 0; mate < quality_by_mate.size() && mate < 2; ++mate) {
        std::array<std::uint32_t, 256> counts{};
        for (const auto base : quality_by_mate[mate])
            ++counts[static_cast<unsigned char>(base)];
        for (std::size_t value = 0; value < counts.size(); ++value)
            if (counts[value] != 0)
                result.by_mate[mate].push_back(QualityCount{
                    static_cast<std::uint8_t>(value), counts[value]});
    }
    return result;
}

std::string formatStarFinalStatistics(const StatsDelta& stats,
    std::chrono::system_clock::time_point started,
    std::chrono::system_clock::time_point mapping_started,
    std::chrono::system_clock::time_point finished) {
    const auto timestamp = [](std::chrono::system_clock::time_point time) {
        const auto raw = std::chrono::system_clock::to_time_t(time);
        const auto local = *std::localtime(&raw);
        std::ostringstream out;
        out << std::put_time(&local, "%b %d %H:%M:%S");
        return out.str();
    };
    const auto fraction = [](double numerator, double denominator) {
        return denominator == 0 ? 0.0 : numerator / denominator;
    };
    const auto percentage = [&](double numerator, double denominator) {
        return fraction(numerator, denominator) * 100.0;
    };
    const auto elapsed = std::difftime(
        std::chrono::system_clock::to_time_t(finished),
        std::chrono::system_clock::to_time_t(mapping_started));
    const auto& splice = stats.splices_by_motif;
    std::ostringstream out;
    out << std::fixed << std::setprecision(2);
    const auto value = [&](const char* label, auto number) {
        out << std::setw(50) << label << number << '\n';
    };
    const auto percent = [&](const char* label, double number) {
        out << std::setw(50) << label << number << "%\n";
    };
    value("Started job on |\t", timestamp(started));
    value("Started mapping on |\t", timestamp(mapping_started));
    value("Finished on |\t", timestamp(finished));
    value("Mapping speed, Million of reads per hour |\t",
          elapsed == 0 ? (stats.reads ? std::numeric_limits<double>::infinity() :
                          std::numeric_limits<double>::quiet_NaN()) :
          double(stats.reads) * 0.0036 / elapsed);
    out << '\n';
    value("Number of input reads |\t", stats.reads);
    value("Average input read length |\t",
          stats.reads ? stats.read_bases / stats.reads : 0);
    out << std::setw(50) << "UNIQUE READS:\n";
    value("Uniquely mapped reads number |\t", stats.mapped_unique);
    percent("Uniquely mapped reads % |\t",
            percentage(stats.mapped_unique, stats.reads));
    value("Average mapped length |\t",
          fraction(stats.mapped_bases, stats.mapped_unique));
    value("Number of splices: Total |\t",
          std::accumulate(splice.begin(), splice.end(), std::uint64_t{0}));
    value("Number of splices: Annotated (sjdb) |\t", stats.annotated_splices);
    value("Number of splices: GT/AG |\t", splice[1] + splice[2]);
    value("Number of splices: GC/AG |\t", splice[3] + splice[4]);
    value("Number of splices: AT/AC |\t", splice[5] + splice[6]);
    value("Number of splices: Non-canonical |\t", splice[0]);
    percent("Mismatch rate per base, % |\t",
            stats.mapped_bases == 0 ?
                std::numeric_limits<double>::quiet_NaN() :
                percentage(stats.mismatches, stats.mapped_bases));
    percent("Deletion rate per base |\t",
            percentage(stats.deletion_bases, stats.mapped_bases));
    value("Deletion average length |\t",
          fraction(stats.deletion_bases, stats.deletion_count));
    percent("Insertion rate per base |\t",
            percentage(stats.insertion_bases, stats.mapped_bases));
    value("Insertion average length |\t",
          fraction(stats.insertion_bases, stats.insertion_count));
    out << std::setw(50) << "MULTI-MAPPING READS:\n";
    value("Number of reads mapped to multiple loci |\t", stats.mapped_multiple);
    percent("% of reads mapped to multiple loci |\t",
            percentage(stats.mapped_multiple, stats.reads));
    value("Number of reads mapped to too many loci |\t", stats.unmapped_multi);
    percent("% of reads mapped to too many loci |\t",
            percentage(stats.unmapped_multi, stats.reads));
    out << std::setw(50) << "UNMAPPED READS:\n";
    value("Number of reads unmapped: too many mismatches |\t", stats.unmapped_mismatch);
    percent("% of reads unmapped: too many mismatches |\t",
            percentage(stats.unmapped_mismatch, stats.reads));
    value("Number of reads unmapped: too short |\t", stats.unmapped_short);
    percent("% of reads unmapped: too short |\t",
            percentage(stats.unmapped_short, stats.reads));
    value("Number of reads unmapped: other |\t", stats.unmapped_other);
    percent("% of reads unmapped: other |\t",
            percentage(stats.unmapped_other, stats.reads));
    out << std::setw(50) << "CHIMERIC READS:\n";
    value("Number of chimeric reads |\t", 0);
    percent("% of chimeric reads |\t", 0);
    return out.str();
}

} // namespace nustar::star2711b
