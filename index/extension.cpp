#include "extension.hpp"

#include <unordered_set>

namespace nustar::star2711b {

IndexExtensionLease SearchIndexBundle::findExtension(std::string_view backend_id) const noexcept {
    for (const auto& extension : extensions)
        if (extension->identity().backend_id == backend_id) return extension;
    return {};
}

Result<SearchIndexBundle> bindIndexExtensions(SearchIndexLease core,
    ReferenceLease reference, std::span<const IndexExtensionLease> extensions) {
    auto invalid = [](std::string message) -> Result<SearchIndexBundle> {
        return Failure{ErrorKind::Index, Stage::Reference, std::nullopt, 1,
                       std::move(message)};
    };
    if (!core || !reference) return invalid("index bundle requires core and reference");
    const auto& core_id = core->identity();
    if (core_id.artifact_sha256.empty() ||
        core_id.reference_semantic_id != reference->view().semantic_id)
        return invalid("core index and reference identity differ");
    SearchIndexBundle bundle{std::move(core), std::move(reference), {}};
    bundle.extensions.reserve(extensions.size());
    std::unordered_set<std::string> names;
    for (const auto& extension : extensions) {
        if (!extension) return invalid("null index extension");
        const auto& identity = extension->identity();
        if (identity.backend_id.empty() || identity.format_id.empty() ||
            identity.format_version == 0 || identity.artifact_sha256.empty() ||
            identity.build_options_sha256.empty())
            return invalid("index extension has incomplete identity");
        if (identity.core_index_sha256 != core_id.artifact_sha256 ||
            identity.reference_semantic_id != core_id.reference_semantic_id)
            return invalid("index extension was built for another core or reference");
        if (!names.insert(identity.backend_id).second)
            return invalid("duplicate index extension backend: " + identity.backend_id);
        bundle.extensions.push_back(extension);
    }
    return bundle;
}

} // namespace nustar::star2711b
