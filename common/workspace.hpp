#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include "types.hpp"

namespace nustar::star2711b {

// Workspace identity includes backend layout and the bound reference/configuration.
struct WorkspaceIdentity {
    std::string backend_id;
    std::uint32_t layout_version;
    std::string binding_fingerprint;
};
struct WorkspaceBudget {
    std::uint64_t scratch_bytes_max;
    std::uint64_t hit_page_size;
};
class StageWorkspace {
public:
    // Execution calls require the workspace created by the same bound kernel
    // and configuration. identity() is metadata, not a per-read runtime check.
    virtual ~StageWorkspace() = default;
    virtual const WorkspaceIdentity& identity() const noexcept = 0;
    virtual void reset() noexcept = 0;
};
class IStageKernel {
public:
    virtual ~IStageKernel() = default;
    virtual std::string_view backendId() const noexcept = 0;
    // No generic factory taking ResolvedRunConfig: each stage declares its actual inputs.
};

} // namespace nustar::star2711b
