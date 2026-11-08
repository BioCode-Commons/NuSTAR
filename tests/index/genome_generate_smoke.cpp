#include "../../runtime/reference_service.hpp"
#include "../../annotation/junction_database.hpp"
#include "../../reference/genome.hpp"
#include "../../reference/synthetic_junction.hpp"
#include "../../index/suffix_array_index.hpp"
#include "../../index/extension.hpp"
#include "../../common/binary.hpp"
#include "../../common/file_identity.hpp"
#include "../../common/sha256.hpp"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

using namespace nustar::star2711b;

struct TestExtension final : IIndexExtension {
    explicit TestExtension(IndexExtensionIdentity identity) : value(std::move(identity)) {}
    const IndexExtensionIdentity& identity() const noexcept override { return value; }
    IndexExtensionIdentity value;
};

struct NaiveSuffixArrayBuilder final : ISuffixArrayBuilder {
    std::string_view algorithmId() const noexcept override { return "naive-test-v1"; }
    std::vector<std::uint32_t> build(
        std::span<const std::uint8_t> symbols) const override {
        std::vector<std::uint32_t> suffixes(symbols.size()+1);
        std::iota(suffixes.begin(), suffixes.end(), 0);
        std::sort(suffixes.begin(), suffixes.end(), [&](auto left, auto right) {
            return std::lexicographical_compare(symbols.begin()+left, symbols.end(),
                                                symbols.begin()+right, symbols.end());
        });
        return suffixes;
    }
};

template<class T> T checked(Result<T> result) {
    if (std::holds_alternative<Failure>(result))
        throw std::runtime_error(std::get<Failure>(result).message);
    return std::get<T>(std::move(result));
}

int main() {
    const auto directory = std::filesystem::temp_directory_path() /
        ("nustar-genome-smoke-" + std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count()));
    struct Cleanup {
        std::filesystem::path path;
        ~Cleanup() { std::error_code error; std::filesystem::remove_all(path, error); }
    } cleanup{directory};
    try {
        Sha256 hash;
        hash.update("abc");
        if (hash.finish() !=
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad")
            throw std::runtime_error("SHA256 implementation failed known vector");
        std::filesystem::create_directories(directory);
        {
            const auto edge_fasta = directory / "padding.fa";
            std::ofstream out(edge_fasta);
            out << ">guard\n" << std::string(200, 'A')
                << "\n>edge\n" << std::string(500, 'C') << '\n';
            out.close();
            const std::array paths{edge_fasta};
            auto storage = checked(readGenomeFasta(paths, 8));
            auto reference = std::make_shared<GenomeReferenceData>(storage,
                std::vector<SjdbEntry>{}, std::vector<BaseCode>{}, 16);
            const auto start = storage->metadata.chromosomes[1].start.value;
            const std::array proposals{
                JunctionProposal{ChromosomeId{1},
                    {GenomeOffset{start + 8}, GenomeOffset{start + 100}},
                    Strand::Plus, JunctionOrigin::JunctionFile, 0, 0, {}},
                JunctionProposal{ChromosomeId{1},
                    {GenomeOffset{start + 350}, GenomeOffset{start + 492}},
                    Strand::Plus, JunctionOrigin::JunctionFile, 0, 1, {}}};
            AnnotationConfig config{};
            config.sjdbOverhang = 16;
            config.limitSjdbInsertNsj = 10;
            GenomeJunctionPreparer preparer;
            auto prepared = checked(preparer.prepare(proposals, reference->view(), config));
            auto sequence = makeSyntheticJunctionSequence(prepared.entries, storage->forward, 16);
            if (prepared.entries.size() != 2 || sequence.size() != 66)
                throw std::runtime_error("chromosome-edge junctions were lost or shortened");
            for (std::size_t i = 0; i < 8; ++i)
                if (sequence[i] != BaseCode::GenomePadding ||
                    sequence[33 + 24 + i] != BaseCode::GenomePadding)
                    throw std::runtime_error("junction flanks did not preserve chromosome padding");
            if (sequence[8] != BaseCode::C || sequence[33 + 23] != BaseCode::C ||
                sequence[32] != BaseCode::GenomePadding || sequence[65] != BaseCode::GenomePadding)
                throw std::runtime_error("junction sequence boundaries changed");
        }
        const auto fasta = directory / "tiny.fa", gtf = directory / "tiny.gtf";
        {
            std::ofstream out(fasta);
            out << ">chr1\n";
            for (int i = 0; i < 500; ++i) out << "ACGT"[i%4];
            out << '\n';
        }
        {
            std::ofstream out(gtf);
            out << "chr1\ttest\texon\t101\t140\t.\t+\t.\tgene_id \"g\"; transcript_id \"t\";\n"
                   "chr1\ttest\texon\t201\t240\t.\t+\t.\tgene_id \"g\"; transcript_id \"t\";\n"
                   "chr1\ttest\texon\t50\t60\t.\t-\t.\tgene_id \"minus\"; transcript_id \"minus\";\n"
                   "chr1\ttest\texon\t50\t60\t.\t+\t.\tgene_id \"plus\"; transcript_id \"plus\";\n";
        }
        GenomeReferenceService references;
        auto base = checked(references.buildBase(ReferenceBuildRequest{{fasta}, directory,
            BackendOptions{"nustar-reference-v1", 1,
                {OptionProvenance{"genomeChrBinNbits", {}, {"8"}, "test"}}}}));
        AnnotationConfig config{};
        config.sjdbGTFfile = gtf;
        config.sjdbGTFfeatureExon = "exon";
        config.sjdbGTFtagExonParentTranscript = "transcript_id";
        config.sjdbGTFtagExonParentGene = "gene_id";
        config.sjdbOverhang = 16;
        config.limitSjdbInsertNsj = 100;
        auto annotated = checked(references.applyAnnotations(base.lease,
            ReferenceAnnotationRequest{config, directory}));
        {
            std::ifstream rows(directory / "exonGeTrInfo.tab");
            std::string count, plus, minus;
            std::getline(rows, count);
            std::getline(rows, plus);
            std::getline(rows, minus);
            if (count != "4" || plus != "49\t59\t1\t2\t2" || minus != "49\t59\t2\t1\t1")
                throw std::runtime_error("coincident exon rows are not ordered by STAR's five columns");
        }
        auto reference = std::dynamic_pointer_cast<const GenomeReferenceData>(annotated.lease);
        if (!reference || reference->junctions().size() != 1 ||
            reference->synthetic().size() != 33)
            throw std::runtime_error("GTF junction did not enter reference");
        const auto genome=reference->starGenome();
        if (genome.data()!=reference->storage().forward.data() ||
            reference->synthetic().data()!=genome.data()+reference->storage().forward.size() ||
            reference->storage().sequence_storage->size()!=genome.size()+40)
            throw std::runtime_error("reference views do not share a single guarded owner");
        const auto copied=*reference;
        if (copied.starGenome().data()!=genome.data() ||
            copied.storageLease()->sequence_storage!=reference->storageLease()->sequence_storage)
            throw std::runtime_error("copying a reference duplicated sequence storage");
        const auto donor = reference->junctions()[0].donor_sequence_begin.value;
        const auto acceptor = reference->junctions()[0].acceptor_sequence_begin.value;
        for (int i = 0; i < 16; ++i) {
            if (reference->synthetic()[i] != reference->storage().forward[donor+i] ||
                reference->synthetic()[i+16] != reference->storage().forward[acceptor+i])
                throw std::runtime_error("inserted splice sequence differs from reference");
        }
        SearchIndexRequest request{"nustar-sa", directory / "index.sa",
            BackendOptions{"nustar-sa-v3", 3,
                {OptionProvenance{"genomeSAindexNbases", {}, {"6"}, "test"}}}};
        SuffixArrayIndexBackend backend;
        auto unsupported_request = request;
        unsupported_request.options.values.push_back(
            OptionProvenance{"genomeSuffixLengthMax", {}, {"80"}, "test"});
        if (!std::holds_alternative<Failure>(
                backend.build(unsupported_request, annotated.lease)))
            throw std::runtime_error("unsupported STAR suffix limit was accepted");
        auto built = checked(backend.build(request, annotated.lease));
        auto loaded_reference = checked(references.loadBase(
            ReferenceRequest{directory, "nustar-reference-v1", {}}));
        auto loaded = checked(backend.load(request, loaded_reference.lease));
        auto index = std::dynamic_pointer_cast<const SuffixArrayIndex>(loaded);
        auto built_index = std::dynamic_pointer_cast<const SuffixArrayIndex>(built);
        if (!index || !built_index || index->identity().format_version != 4 ||
            index->suffixArray().bitWidth() >= 32 ||
            index->suffixArray().packedBytes().size() >=
                index->suffixArray().size()*sizeof(std::uint32_t) ||
            !std::ranges::equal(index->suffixArray(), built_index->suffixArray()))
            throw std::runtime_error("packed SA reload or compression differs");
        for(const auto algorithm:{"libsais64","libsais40"}) {
            auto selected=request;selected.artifact=directory/(std::string(algorithm)+".sa");
            selected.options.values.push_back(OptionProvenance{"genomeSAalgorithm",{},{algorithm},"test"});
            auto other=std::dynamic_pointer_cast<const SuffixArrayIndex>(checked(backend.build(selected,annotated.lease)));
            auto cross=selected;cross.options.values.back().resolved_tokens={"libsais64"};
            auto reloaded=std::dynamic_pointer_cast<const SuffixArrayIndex>(checked(backend.load(cross,loaded_reference.lease)));
            if(!std::ranges::equal(index->suffixArray(),other->suffixArray()) ||
               !std::ranges::equal(index->suffixArray(),reloaded->suffixArray()) ||
               !std::ranges::equal(index->prefixIndex().bytes(),other->prefixIndex().bytes()))
                throw std::runtime_error("selected SA backend changed index content");
        }
        if (backend.builderId()!="libsais64-2.12.1")
            throw std::runtime_error("default builder provenance differs");
        for (const auto algorithm:{"unknown","sais32","psacak40"}) {
            auto unknown=request;
            unknown.options.values.push_back(OptionProvenance{"genomeSAalgorithm",{},{algorithm},"test"});
            if(!std::holds_alternative<Failure>(backend.build(unknown,annotated.lease)))
                throw std::runtime_error("unsupported SA backend accepted");
        }
        // The sidecar is optional build provenance, never a load-time check.
        {
            std::ofstream bad_checksum(request.artifact.string()+".sha256");
            bad_checksum << std::string(64,'0') << '\n';
        }
        checked(backend.load(request, loaded_reference.lease));
        std::filesystem::remove(request.artifact.string()+".sha256");
        auto without_sidecar=std::dynamic_pointer_cast<const SuffixArrayIndex>(
            checked(backend.load(request, loaded_reference.lease)));
        if (!std::ranges::equal(index->suffixArray(),without_sidecar->suffixArray()))
            throw std::runtime_error("optional provenance changed loaded SA");
        {
            std::ofstream good_checksum(request.artifact.string()+".sha256");
            good_checksum << index->identity().artifact_sha256 << '\n';
        }
        const auto splice_position = 2 * reference->storage().forward.size() + 40;
        // SA entries are ordered by suffix content, not by their numeric position.
        if (!index || std::ranges::find(index->suffixArray(), splice_position) ==
            index->suffixArray().end())
            throw std::runtime_error("splice sequence missing from SAIS index");
        // A second algorithm must produce the same STAR-normalized index and remain
        // loadable through the default backend without knowing that algorithm.
        auto alternate_request = request;
        alternate_request.artifact = directory / "index.naive-test.sa";
        SuffixArrayIndexBackend alternate(std::make_shared<NaiveSuffixArrayBuilder>());
        checked(alternate.build(alternate_request, annotated.lease));
        auto alternate_loaded = std::dynamic_pointer_cast<const SuffixArrayIndex>(
            checked(backend.load(alternate_request, loaded_reference.lease)));
        if (!alternate_loaded || alternate_loaded->algorithmId() != "naive-test-v1" ||
            !std::ranges::equal(index->suffixArray(), alternate_loaded->suffixArray()) ||
            !std::ranges::equal(index->prefixRanges(), alternate_loaded->prefixRanges(),
                [](const PrefixRange& left, const PrefixRange& right) {
                    return left.code == right.code && left.begin == right.begin &&
                           left.end == right.end;
                }))
            throw std::runtime_error("suffix array builder substitution changed the index");
        // A v2 artifact made from the same ranks must still load through the v3
        // backend. This keeps existing user indexes usable after the upgrade.
        auto legacy_request = request;
        legacy_request.artifact = directory / "index.legacy-v2.sa";
        legacy_request.options.schema_id = "nustar-sa-v2";
        legacy_request.options.schema_version = 2;
        Sha256 legacy_options;
        legacy_options.update("nustar-sa-v2:");
        legacy_options.update("6:");
        legacy_options.update(std::to_string(std::numeric_limits<std::uint32_t>::max()));
        legacy_options.update(":");
        legacy_options.update(index->algorithmId());
        {
            std::ofstream out(legacy_request.artifact, std::ios::binary);
            out.write("KSAI", 4);
            binary::put32(out, 2);
            binary::put32(out, index->prefixBases());
            binary::put32(out, index->suffixLengthMax());
            binary::put64(out, index->physicalLength());
            binary::put64(out, index->suffixArray().size());
            binary::put64(out, index->prefixRanges().size());
            binary::putString(out, index->identity().reference_semantic_id);
            binary::putString(out, legacy_options.finish());
            binary::putString(out, std::string(index->algorithmId()));
            for (const auto position : index->suffixArray())
                binary::put32(out, position);
            for (const auto& prefix : index->prefixRanges()) {
                binary::put64(out, prefix.code);
                binary::put64(out, prefix.begin);
                binary::put64(out, prefix.end);
            }
            if (!out) throw std::runtime_error("cannot write legacy index fixture");
        }
        const auto legacy_file = checked(identifyFile(legacy_request.artifact));
        {
            std::ofstream checksum(legacy_request.artifact.string()+".sha256");
            checksum << legacy_file.sha256 << '\n';
        }
        auto legacy_loaded = std::dynamic_pointer_cast<const SuffixArrayIndex>(
            checked(backend.load(legacy_request, loaded_reference.lease)));
        if (!legacy_loaded || legacy_loaded->identity().format_version != 2 ||
            !std::ranges::equal(index->suffixArray(), legacy_loaded->suffixArray()))
            throw std::runtime_error("legacy v2 suffix array did not reload");
        auto extension = std::make_shared<TestExtension>(IndexExtensionIdentity{
            "lcp-test", "lcp-v1", 1, "extension-hash",
            loaded->identity().artifact_sha256,
            loaded->identity().reference_semantic_id, "option-hash"});
        const std::array<IndexExtensionLease, 1> matching{extension};
        auto bundle = checked(bindIndexExtensions(loaded, loaded_reference.lease, matching));
        if (bundle.findExtension("lcp-test") != extension)
            throw std::runtime_error("valid index extension did not bind");
        const std::array<IndexExtensionLease, 1> stale{
            std::make_shared<TestExtension>(IndexExtensionIdentity{
                "lcp-test", "lcp-v1", 1, "extension-hash", "stale-core",
                loaded->identity().reference_semantic_id, "option-hash"})};
        if (!std::holds_alternative<Failure>(
                bindIndexExtensions(loaded, loaded_reference.lease, stale)))
            throw std::runtime_error("stale index extension was accepted");
        std::fstream artifact(request.artifact, std::ios::in | std::ios::out | std::ios::binary);
        artifact.seekp(100); artifact.put(char(255)); artifact.close();
        if (!std::holds_alternative<Failure>(backend.load(request, loaded_reference.lease)))
            throw std::runtime_error("corrupted index was accepted");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
