#pragma once

#include "../reference/genome.hpp"
#include "../annotation/config.hpp"
#include "../common/backend_options.hpp"

namespace nustar::star2711b {

struct ReferenceRequest {
    std::filesystem::path source; // A STAR bundle or another declared reference format.
    std::string source_format;
    BackendOptions storage_options; // Load the bundle's existing annotation state only.
};
struct ReferenceBuildRequest {
    std::vector<std::filesystem::path> fasta_files;
    std::filesystem::path output_directory;
    BackendOptions storage_options; // Reference storage only; search-index options are separate.
};
struct ReferenceAnnotationRequest {
    AnnotationConfig annotation; // Only runtime or genomeGenerate additions to the base reference.
    std::filesystem::path output_directory;
};
struct ReferenceReady {
    ReferenceLease lease;
    std::vector<FileIdentity> files;
    std::vector<Diagnostic> diagnostics;
};
class IReferenceService {
public:
    virtual ~IReferenceService() = default;
    // Existing index annotations belong to this immutable base generation.
    virtual Result<ReferenceReady> loadBase(const ReferenceRequest&) = 0;
    // Build sequence/metadata from FASTA without runtime annotation additions.
    virtual Result<ReferenceReady> buildBase(const ReferenceBuildRequest&) = 0;
    // Publish a new immutable generation; the base remains alive for old index readers.
    virtual Result<ReferenceReady> applyAnnotations(
        ReferenceLease base, const ReferenceAnnotationRequest&) = 0;
};

class GenomeReferenceService final : public IReferenceService {
public:
    Result<ReferenceReady> loadBase(const ReferenceRequest&) override;
    Result<ReferenceReady> buildBase(const ReferenceBuildRequest&) override;
    Result<ReferenceReady> applyAnnotations(
        ReferenceLease base, const ReferenceAnnotationRequest&) override;
};

} // namespace nustar::star2711b
