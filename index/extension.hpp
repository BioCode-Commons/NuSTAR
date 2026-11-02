#pragma once

#include "backend.hpp"

namespace nustar::star2711b {

// Auxiliary search structures (for example LCP or NuSTAR's extended k-mer hash)
// are separate artifacts. Rebuilding one never changes the core SA format.
struct IndexExtensionIdentity {
    std::string backend_id, format_id;
    std::uint32_t format_version;
    std::string artifact_sha256, core_index_sha256, reference_semantic_id,
                build_options_sha256;
};

class IIndexExtension {
public:
    virtual ~IIndexExtension() = default;
    virtual const IndexExtensionIdentity& identity() const noexcept = 0;
};
using IndexExtensionLease = std::shared_ptr<const IIndexExtension>;

struct IndexExtensionRequest {
    std::string backend_id;
    std::filesystem::path artifact;
    BackendOptions options;
};

class IIndexExtensionBackend {
public:
    virtual ~IIndexExtensionBackend() = default;
    virtual std::string_view backendId() const noexcept = 0;
    // The concrete extension can downcast a compatible core index and use its
    // native SA/prefix views. The common ISearchIndex remains format agnostic.
    virtual Result<IndexExtensionLease> build(const IndexExtensionRequest&,
        SearchIndexLease core, ReferenceLease reference) = 0;
    virtual Result<IndexExtensionLease> load(const IndexExtensionRequest&,
        SearchIndexLease core, ReferenceLease reference) = 0;
    virtual Status save(const IIndexExtension&, const std::filesystem::path&) = 0;
};

struct SearchIndexBundle {
    SearchIndexLease core;
    ReferenceLease reference;
    std::vector<IndexExtensionLease> extensions;
    IndexExtensionLease findExtension(std::string_view backend_id) const noexcept;
};

// Publish the bundle only when every extension names this exact core artifact
// and reference generation. This prevents stale LCP/hash tables after GTF updates.
Result<SearchIndexBundle> bindIndexExtensions(SearchIndexLease core,
    ReferenceLease reference, std::span<const IndexExtensionLease> extensions);

} // namespace nustar::star2711b
