#pragma once

#include "../common/workspace.hpp"
#include "../pipeline/effects.hpp"
#include "../reads/model.hpp"
#include "../reference/view.hpp"
#include "../selection/model.hpp"
#include "config.hpp"
#include <filesystem>
#include <memory>

namespace nustar::star2711b::solo {

Status filterCells(const Config&, const std::filesystem::path& input,
                   const std::string& output_prefix, const std::string& log_prefix);

class Worker {
  public:
    struct Impl;
    explicit Worker(std::shared_ptr<Impl>);
    void record(const RawMate* barcode, const MappingRead&, const CandidateSet&,
                const SelectionResult&, const FilterDecision&, const StatsDelta&,
                const QualityHistogramDelta&, ProjectedOutput* = nullptr,
                const std::string& input_tags = {});

  private:
    std::shared_ptr<Impl> impl_;
};

class Session {
  public:
    struct Impl;
    Session(const Config&, const std::filesystem::path& genome_dir, std::string prefix,
            unsigned threads, unsigned mates, const ReferenceMetadata&);
    ~Session();
    std::unique_ptr<Worker> createWorker();
    void finish();
    bool hasCorrections() const;
    std::pair<std::string, std::string> correctedTags(std::uint64_t read) const;

  private:
    std::unique_ptr<Impl> impl_;
};

// Runtime owns this composition; the candidate mapper has no Solo dependency.
class Workspace final : public StageWorkspace {
  public:
    Workspace(std::unique_ptr<StageWorkspace> mapping, std::unique_ptr<Worker> solo)
        : mapping(std::move(mapping)), counting(std::move(solo)) {}
    const WorkspaceIdentity& identity() const noexcept override { return mapping->identity(); }
    void reset() noexcept override { mapping->reset(); }
    std::unique_ptr<StageWorkspace> mapping;
    std::unique_ptr<Worker> counting;
};

} // namespace nustar::star2711b::solo
