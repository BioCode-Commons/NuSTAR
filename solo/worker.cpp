#include "worker_state.hpp"
#include <stdexcept>

namespace nustar::star2711b::solo {
using detail::EX_G;
using detail::EX_iFrag;
using detail::EX_L;
using detail::EX_R;
namespace FeatureTypes = detail::FeatureTypes;
using detail::UnsignedWord;
Worker::Worker(std::shared_ptr<Impl> impl) : impl_(std::move(impl)) {}
void Worker::record(const RawMate* barcode, const MappingRead& read, const CandidateSet& candidates,
                    const SelectionResult& selection, const FilterDecision& decision,
                    const StatsDelta& statistics, const QualityHistogramDelta& quality,
                    ProjectedOutput* output, const std::string& input_tags) {
    auto& state = *impl_;
    auto& p = state.parameters;
    addStatistics(state.statistics, statistics);
    for (unsigned m = 0; m < read.mates.size(); ++m)
        for (auto q : quality.by_mate[m])
            state.native.qualHist[m][q.value] += q.count;
    char* sequences[3]{};
    char* qualities[3]{};
    UnsignedWord lengths[3]{};
    if (barcode) {
        sequences[p.options.barcodeRead] = const_cast<char*>(barcode->sequence.data());
        qualities[p.options.barcodeRead] = const_cast<char*>(barcode->quality->data());
        lengths[p.options.barcodeRead] = barcode->sequence.size();
    }
    if (p.readFilesTypeN == 10 && input_tags.empty())
        throw std::runtime_error("Solo SAM/BAM input has no barcode tags");
    auto& solo = *state.native.read;
    solo.readBar->getCBandUMI(sequences, qualities, lengths, input_tags, read.fragment.file_group,
                              "NuSTAR read");
    const auto count = decision.mapping == MappingClass::Unmapped ? 0 : selection.ordered.size();
    state.transcripts.resize(count);
    state.alignments.resize(count);
    for (std::size_t i = 0; i < count; ++i) {
        const auto& selected = selection.ordered[i];
        const auto& c = candidates.candidates[selected.candidate.value];
        auto& t = state.transcripts[i];
        t.nExons = c.blocks.size();
        t.Str = static_cast<unsigned>(c.orientation);
        t.Chr = c.chromosome.value;
        t.cStart = selected.chromosome_relative_start;
        t.Lread = read.budget.legacy_Lread;
        t.sjYes = c.metrics.has_splice;
        t.readNmates = read.mates.size();
        for (std::size_t k = 0; k < c.blocks.size(); ++k) {
            const auto& b = c.blocks[k];
            t.exons[k][EX_R] = b.read_start.value;
            t.exons[k][EX_G] = b.genome_start.value;
            t.exons[k][EX_L] = b.length;
            t.exons[k][EX_iFrag] = static_cast<unsigned>(b.mate);
        }
        for (std::size_t k = 0; k < c.boundaries.size(); ++k) {
            const auto& b = c.boundaries[k];
            t.canonSJ[k] = b.kind == BoundaryKind::Intron      ? static_cast<int>(b.motif)
                           : b.kind == BoundaryKind::Deletion  ? -1
                           : b.kind == BoundaryKind::Insertion ? -2
                                                               : -3;
            t.sjAnnot[k] = b.annotation;
        }
        state.alignments[i] = &t;
    }
    auto& a = state.annotations;
    a.reset();
    auto** align = state.alignments.data();
    auto& tr = state.transcriptome;
    if (p.quant.geneFull.yes)
        tr.geneFullAlignOverlap(count, align, p.options.strand,
                                a.annotFeatures[FeatureTypes::GeneFull]);
    if (p.quant.gene.yes)
        tr.classifyAlign(align, count, a);
    if (p.quant.geneFull_ExonOverIntron.yes)
        tr.geneFullAlignOverlap_ExonOverIntron(
            count, align, p.options.strand, a.annotFeatures[FeatureTypes::GeneFull_ExonOverIntron],
            a.annotFeatures[FeatureTypes::Gene]);
    if (p.quant.geneFull_Ex50pAS.yes)
        tr.alignExonOverlap(count, align, p.options.strand,
                            a.annotFeatures[FeatureTypes::GeneFull_Ex50pAS]);
    solo.record(count, align, read.fragment.ordinal, a);
    if (!output || p.outputAttributes.empty())
        return;
    const auto& bar = *solo.readBar;
    for (auto& record : output->records) {
        auto append = [&](const std::string& tag, char type, std::string value) {
            record.tags.push_back({{tag[0], tag[1]}, type, std::move(value)});
        };
        const auto& feature = a.annotFeatures[p.options.samAttrFeature];
        const auto alignment = record.alignment_ordinal;
        const bool mapped = !(record.flag & 4);
        for (const auto& tag : p.outputAttributes) {
            if (tag == "CR")
                append(tag, 'Z', bar.cbSeq);
            else if (tag == "CY")
                append(tag, 'Z', bar.cbQual);
            else if (tag == "UR")
                append(tag, 'Z', bar.umiSeq);
            else if (tag == "UY")
                append(tag, 'Z', bar.umiQual);
            else if (tag == "CB" && !bar.cbSeqCorrected.empty())
                append(tag, 'Z', bar.cbSeqCorrected);
            else if (tag == "sM")
                append(tag, 'i', std::to_string(bar.cbMatch));
            else if (tag == "sS")
                append(tag, 'Z', bar.bSeq);
            else if (tag == "sQ")
                append(tag, 'Z', bar.bQual);
            else if (mapped && (tag == "GX" || tag == "GN")) {
                std::string value = "-";
                if (feature.fSet.size() == 1 && alignment < feature.fAlign.size() &&
                    feature.fAlign[alignment].size() == 1) {
                    const auto gene = *feature.fAlign[alignment].begin();
                    value = tag == "GX" ? tr.geID[gene] : tr.geName[gene];
                }
                append(tag, 'Z', std::move(value));
            } else if (mapped && tag == "sF") {
                int overlap = feature.ovType, genes = feature.fSet.size();
                if ((overlap == 1 || overlap == 3 || overlap == 5) &&
                    alignment < feature.fAlign.size() && feature.fAlign[alignment].empty())
                    overlap = genes = -1;
                append(tag, 'B', "i," + std::to_string(overlap) + "," + std::to_string(genes));
            }
        }
        {
            static const std::vector<std::string> unmapped_order{"NH", "HI", "AS", "nM", "uT", "RG", "CR",
                                                          "CY", "CB", "UR", "UY", "sM", "sS", "sQ"};
            const auto& order = mapped ? p.outputAttributes : unmapped_order;
            std::vector<SamAuxField> ordered;
            ordered.reserve(record.tags.size());
            for (const auto& tag : order)
                for (auto& value : record.tags)
                    if (value.tag[0] == tag[0] && value.tag[1] == tag[1]) {
                        ordered.push_back(std::move(value));
                        break;
                    }
            record.tags = std::move(ordered);
        }
    }
}

} // namespace nustar::star2711b::solo
