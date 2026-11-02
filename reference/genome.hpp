#pragma once

#include <filesystem>
#include <memory>
#include <span>
#include <string>
#include <vector>
#include "view.hpp"

namespace nustar::star2711b {

struct GenomeStorage {
    ReferenceMetadata metadata;
    // Immutable sequence view and its owner travel together, including copies.
    std::shared_ptr<const std::vector<BaseCode>> sequence_storage;
    std::span<const BaseCode> forward; // Chromosomes plus bin padding.
    std::uint32_t chromosome_bin_bits = 18;
};

class GenomeReferenceData final : public IReferenceData, public IReferenceSequence,
                                  public ISequenceLease,
                                  public std::enable_shared_from_this<GenomeReferenceData> {
public:
    GenomeReferenceData(std::shared_ptr<const GenomeStorage> storage,
        std::vector<SjdbEntry> junctions, std::vector<BaseCode> synthetic,
        std::uint64_t sjdb_overhang, bool transient = false);

    ReferenceView view() const override;
    Result<SequenceWindow> pin(GenomeRange range) const override;

    const GenomeStorage& storage() const noexcept { return *storage_; }
    std::shared_ptr<const GenomeStorage> storageLease() const noexcept { return storage_; }
    std::span<const BaseCode> synthetic() const noexcept { return synthetic_; }
    std::span<const SjdbEntry> junctions() const noexcept { return junctions_; }
    std::uint64_t sjdbOverhang() const noexcept { return sjdb_overhang_; }
    std::string_view semanticId() const noexcept { return semantic_id_; }
    // STAR Genome.G layout: forward chromosomes followed by synthetic SJ bases.
    // One immutable, reference-owned buffer shared by seed and stitch modules.
    // The span remains valid for the lifetime of this reference's lease.
    std::span<const BaseCode> starGenome() const;

private:
    std::shared_ptr<const GenomeStorage> storage_;
    std::vector<SjdbEntry> junctions_;
    std::span<const BaseCode> synthetic_;
    std::uint64_t sjdb_overhang_;
    ReferenceGeneration generation_;
    std::string semantic_id_;
    std::shared_ptr<const std::vector<BaseCode>> star_genome_;
    void identify();
    // Loader already owns one guarded forward+SJ buffer; do not copy it.
    GenomeReferenceData(std::shared_ptr<const GenomeStorage>,
        std::vector<SjdbEntry>, std::uint64_t sjdb_overhang,
        std::shared_ptr<const std::vector<BaseCode>> genome, std::string semantic_id);
    friend Result<std::shared_ptr<const GenomeReferenceData>> loadGenomeReference(
        const std::filesystem::path&);
};

Result<std::shared_ptr<const GenomeStorage>> readGenomeFasta(
    std::span<const std::filesystem::path> files, std::uint32_t chromosome_bin_bits);
Result<FileIdentity> saveGenomeReference(const GenomeReferenceData& reference,
    const std::filesystem::path& directory);
Result<std::shared_ptr<const GenomeReferenceData>> loadGenomeReference(
    const std::filesystem::path& directory);

} // namespace nustar::star2711b
