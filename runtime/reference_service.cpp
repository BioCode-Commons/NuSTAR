#include "reference_service.hpp"
#include "../annotation/table_writer.hpp"
#include "../annotation/junction_database.hpp"
#include "../reference/synthetic_junction.hpp"
#include "../annotation/gtf_loader.hpp"
#include "../common/file_identity.hpp"

#include <charconv>

namespace nustar::star2711b {
namespace {

Failure serviceError(ErrorKind kind, Stage stage, std::string message) {
    return Failure{kind, stage, std::nullopt, 1, std::move(message)};
}

Result<std::uint32_t> binBits(const BackendOptions& options) {
    if (!options.schema_id.empty() &&
        (options.schema_id != "nustar-reference-v1" || options.schema_version != 1))
        return serviceError(ErrorKind::Config, Stage::Reference, "unsupported reference options schema");
    std::uint32_t bits = 18;
    for (const auto& value : options.values) {
        if (value.name != "genomeChrBinNbits" || value.resolved_tokens.size() != 1)
            return serviceError(ErrorKind::Config, Stage::Reference,
                                "unknown reference option: " + value.name);
        const auto& token = value.resolved_tokens[0];
        std::uint32_t parsed = 0;
        const auto [end, error] = std::from_chars(token.data(), token.data()+token.size(), parsed);
        if (error != std::errc{} || end != token.data()+token.size() ||
            parsed < 8 || parsed > 30) {
            return serviceError(ErrorKind::Config, Stage::Reference,
                                "invalid genomeChrBinNbits");
        }
        bits = parsed;
    }
    return bits;
}

} // namespace

Result<ReferenceReady> GenomeReferenceService::buildBase(const ReferenceBuildRequest& request) {
    auto bits = binBits(request.storage_options);
    if (std::holds_alternative<Failure>(bits)) return std::get<Failure>(std::move(bits));
    auto storage = readGenomeFasta(request.fasta_files, std::get<std::uint32_t>(bits));
    if (std::holds_alternative<Failure>(storage)) return std::get<Failure>(std::move(storage));
    auto reference = std::make_shared<const GenomeReferenceData>(
        std::get<std::shared_ptr<const GenomeStorage>>(std::move(storage)),
        std::vector<SjdbEntry>{}, std::vector<BaseCode>{}, 0);
    auto saved = saveGenomeReference(*reference, request.output_directory);
    if (std::holds_alternative<Failure>(saved)) return std::get<Failure>(std::move(saved));
    return ReferenceReady{reference, {std::get<FileIdentity>(std::move(saved))}, {}};
}

Result<ReferenceReady> GenomeReferenceService::loadBase(const ReferenceRequest& request) {
    if (request.source_format != "nustar-reference-v1")
        return serviceError(ErrorKind::Config, Stage::Reference,
                            "unsupported reference format: " + request.source_format);
    auto loaded = loadGenomeReference(request.source);
    if (std::holds_alternative<Failure>(loaded)) return std::get<Failure>(std::move(loaded));
    auto identity = identifyFile(request.source / "reference.bin");
    if (std::holds_alternative<Failure>(identity)) return std::get<Failure>(std::move(identity));
    return ReferenceReady{
        std::get<std::shared_ptr<const GenomeReferenceData>>(std::move(loaded)),
        {std::get<FileIdentity>(std::move(identity))}, {}};
}

Result<ReferenceReady> GenomeReferenceService::applyAnnotations(
    ReferenceLease base, const ReferenceAnnotationRequest& request) {
    const auto concrete = std::dynamic_pointer_cast<const GenomeReferenceData>(base);
    if (!concrete) return serviceError(ErrorKind::Config, Stage::Reference,
                                       "annotations require a NuSTAR reference");
    if (!request.annotation.sjdbFileChrStartEnd.empty())
        return serviceError(ErrorKind::Config, Stage::Sjdb,
                            "sjdbFileChrStartEnd is not implemented for this backend");
    if (!request.annotation.sjdbGTFfile)
        return serviceError(ErrorKind::Config, Stage::Gtf, "no GTF annotation was supplied");
    GtfAnnotationLoader loader;
    auto catalog_result = loader.load(request.annotation, concrete->view());
    if (std::holds_alternative<Failure>(catalog_result))
        return std::get<Failure>(std::move(catalog_result));
    auto catalog = std::get<AnnotationCatalog>(std::move(catalog_result));
    auto proposal_result = loader.deriveJunctions(catalog, concrete->view());
    if (std::holds_alternative<Failure>(proposal_result))
        return std::get<Failure>(std::move(proposal_result));
    auto proposals = std::get<std::vector<JunctionProposal>>(std::move(proposal_result));
    GtfTableWriter writer;
    auto status = writer.writeLegacyTables(catalog, proposals,
        concrete->storage().metadata, request.output_directory);
    if (std::holds_alternative<Failure>(status)) return std::get<Failure>(std::move(status));
    GenomeJunctionPreparer preparer;
    auto prepared_result = preparer.prepare(proposals, concrete->view(), request.annotation);
    if (std::holds_alternative<Failure>(prepared_result))
        return std::get<Failure>(std::move(prepared_result));
    auto prepared = std::get<PreparedJunctionDatabase>(std::move(prepared_result));
    const auto overhang = prepared.overhang;
    auto synthetic = makeSyntheticJunctionSequence(prepared.entries,
        concrete->storage().forward, overhang);
    auto annotated = std::make_shared<const GenomeReferenceData>(concrete->storageLease(),
        std::move(prepared.entries), std::move(synthetic), overhang);
    auto saved = saveGenomeReference(*annotated, request.output_directory);
    if (std::holds_alternative<Failure>(saved)) return std::get<Failure>(std::move(saved));
    status = writeJunctionDatabaseTables(annotated->view(), request.output_directory);
    if (std::holds_alternative<Failure>(status)) return std::get<Failure>(std::move(status));
    std::vector<FileIdentity> files{std::get<FileIdentity>(std::move(saved))};
    for (const auto& filename : {"geneInfo.tab", "exonGeTrInfo.tab", "transcriptInfo.tab",
                                 "exonInfo.tab", "sjdbList.fromGTF.out.tab",
                                 "sjdbInfo.txt", "sjdbList.out.tab"}) {
        auto identity = identifyFile(request.output_directory / filename);
        if (std::holds_alternative<Failure>(identity)) return std::get<Failure>(std::move(identity));
        files.push_back(std::get<FileIdentity>(std::move(identity)));
    }
    auto diagnostics = std::move(catalog.diagnostics);
    diagnostics.insert(diagnostics.end(), prepared.diagnostics.begin(), prepared.diagnostics.end());
    return ReferenceReady{annotated, std::move(files), std::move(diagnostics)};
}

} // namespace nustar::star2711b
