#pragma once

#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <string_view>
#include <vector>
#include "../common/backend_options.hpp"
#include "../reference/view.hpp"

namespace nustar::star2711b {

// Search indexes and seed algorithms are selected independently, then explicitly bound.
struct SearchIndexIdentity {
    std::string backend_id, format_id;
    std::uint32_t format_version;
    // Artifact digest is optional build provenance, not verified during loading.
    std::string artifact_sha256, reference_semantic_id, build_options_sha256;
    IndexGeneration generation;
};
struct SearchIndexCapabilities {
    std::vector<std::string> native_query_contracts; // Versioned interfaces understood by binders.
    bool incremental_annotation_update;
};
class ISearchIndex {
public:
    virtual ~ISearchIndex() = default;
    virtual const SearchIndexIdentity& identity() const noexcept = 0;
    virtual SearchIndexCapabilities capabilities() const = 0;
    // No mandatory SA, rank, LF, minimizer, or longest-prefix method belongs here.
};
using SearchIndexLease = std::shared_ptr<const ISearchIndex>;
struct SearchIndexRequest {
    std::string backend_id;
    std::filesystem::path artifact;
    BackendOptions options;
};
class IIndexBackend {
public:
    virtual ~IIndexBackend() = default;
    virtual std::string_view backendId() const noexcept = 0;
    // Load against the base reference encoded by the artifact, before runtime annotations.
    virtual Result<SearchIndexLease> load(const SearchIndexRequest&, ReferenceLease base) = 0;
    virtual Result<SearchIndexLease> build(const SearchIndexRequest&, ReferenceLease) = 0;
    virtual Status save(const ISearchIndex&, const std::filesystem::path&) = 0;
    // previous retains its base reference; updated_reference has a distinct semantic identity.
    // Explicitly fail if unsupported; the planner may choose a full rebuild before mapping.
    virtual Result<SearchIndexLease> updateAnnotations(
        SearchIndexLease previous, ReferenceLease updated_reference,
        const BackendOptions&) = 0;
};

} // namespace nustar::star2711b
