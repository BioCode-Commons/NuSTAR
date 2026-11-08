#pragma once

#include <utility>
#include <cstdint>
#include <filesystem>
#include <span>
#include <string>
#include <vector>
#include "suffix_array.hpp"
#include "backend.hpp"
#include "packed_suffix_array.hpp"

namespace nustar::star2711b {

class GenomeReferenceData;

struct PrefixRange {
    std::uint64_t code, begin, end; // Half-open SA range for a full A/C/G/T prefix.
};

// STAR multilevel SAi: direct-address ranks plus N/absent flags. The legacy
// sparse table is accepted only at the disk boundary, never kept during mapping.
class StarPrefixIndex {
public:
    StarPrefixIndex() = default;
    StarPrefixIndex(std::uint32_t bases, std::uint64_t genome_length,
        const PackedSuffixArray&, const GenomeReferenceData&,
        std::span<const PrefixRange> legacy);
    static StarPrefixIndex fromBytes(std::uint32_t bases, std::uint64_t genome_length,
        std::vector<std::uint8_t> bytes);
    static StarPrefixIndex fromView(std::uint32_t bases, std::uint64_t genome_length,
        std::span<const std::uint8_t> bytes,
        std::shared_ptr<const void> owner);
    static std::uint64_t entryCount(std::uint32_t bases);
    static unsigned width(std::uint64_t genome_length);
    std::uint64_t operator[](std::uint64_t i) const noexcept { return entries_[i]; }
    std::uint64_t levelStart(unsigned length) const noexcept { return starts_[length-1]; }
    std::uint64_t levelEnd(unsigned length) const noexcept { return starts_[length]; }
    std::uint64_t nMask() const noexcept { return n_mask_; }
    std::uint64_t absentMask() const noexcept { return absent_mask_; }
    std::uint64_t rankMask() const noexcept { return n_mask_-1; }
    std::span<const std::uint8_t> bytes() const noexcept { return entries_.packedBytes(); }
private:
    void layout(unsigned bases, std::uint64_t genome_length);
    std::vector<std::uint64_t> starts_;
    PackedSuffixArray entries_;
    std::uint64_t n_mask_ = 0, absent_mask_ = 0;
};

class SuffixArrayIndex final : public ISearchIndex {
public:
    SuffixArrayIndex(SearchIndexIdentity identity, ReferenceLease reference,
        std::filesystem::path artifact, std::string algorithm_id,
        std::uint32_t prefix_bases, std::uint32_t suffix_length_max,
        std::uint64_t physical_length, PackedSuffixArray suffixes,
        std::vector<PrefixRange> prefixes, StarPrefixIndex native_prefix = {});
    const SearchIndexIdentity& identity() const noexcept override { return identity_; }
    SearchIndexCapabilities capabilities() const override;
    const PackedSuffixArray& suffixArray() const noexcept { return suffixes_; }
    // Diagnostic compatibility adapter; materializes the old ranges on request.
    std::vector<PrefixRange> prefixRanges() const;
    const StarPrefixIndex& prefixIndex() const noexcept { return prefix_index_; }
    std::string_view algorithmId() const noexcept { return algorithm_id_; }
    std::uint32_t prefixBases() const noexcept { return prefix_bases_; }
    std::uint32_t suffixLengthMax() const noexcept { return suffix_length_max_; }
    std::uint64_t physicalLength() const noexcept { return physical_length_; }
    const std::filesystem::path& artifact() const noexcept { return artifact_; }
private:
    friend class SuffixArrayIndexBackend;
    void setArtifactSha256(std::string hash) { identity_.artifact_sha256 = std::move(hash); }
    SearchIndexIdentity identity_;
    ReferenceLease reference_;
    std::filesystem::path artifact_;
    std::string algorithm_id_;
    std::uint32_t prefix_bases_;
    std::uint32_t suffix_length_max_;
    std::uint64_t physical_length_;
    PackedSuffixArray suffixes_;
    StarPrefixIndex prefix_index_;
};

class SuffixArrayIndexBackend final : public IIndexBackend {
public:
    explicit SuffixArrayIndexBackend(SuffixArrayBuilderLease builder = {});
    std::string_view backendId() const noexcept override { return "nustar-sa"; }
    std::string_view builderId() const noexcept {
        return builder_ ? builder_->algorithmId() : "libsais64-2.12.1";
    }
    Result<SearchIndexLease> load(const SearchIndexRequest&, ReferenceLease base) override;
    Result<SearchIndexLease> build(const SearchIndexRequest&, ReferenceLease) override;
    Status save(const ISearchIndex&, const std::filesystem::path&) override;
    Result<SearchIndexLease> updateAnnotations(
        SearchIndexLease previous, ReferenceLease updated_reference,
        const BackendOptions&) override;
private:
    SuffixArrayBuilderLease builder_;
};

} // namespace nustar::star2711b
