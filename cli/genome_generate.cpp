#include "../runtime/reference_service.hpp"
#include "genome_generate.hpp"
#include "../reference/genome.hpp"
#include "../index/suffix_array_index.hpp"

#include <charconv>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {
using namespace nustar::star2711b;

std::uint64_t unsignedOption(std::string_view token, std::string_view name,
                             std::uint64_t low, std::uint64_t high) {
    std::uint64_t value = 0;
    const auto result = std::from_chars(token.data(), token.data() + token.size(), value);
    if (result.ec != std::errc{} || result.ptr != token.data() + token.size() ||
        value < low || value > high)
        throw std::invalid_argument("invalid value for " + std::string(name));
    return value;
}

template<class T> T take(Result<T> result) {
    if (std::holds_alternative<Failure>(result))
        throw std::runtime_error(std::get<Failure>(result).message);
    return std::get<T>(std::move(result));
}

} // namespace

int runGenomeGenerate(int argc, char** argv, int first_option) {
    try {
        ReferenceBuildRequest build;
        build.output_directory = "./GenomeDir/";
        AnnotationConfig annotation{};
        annotation.sjdbGTFfeatureExon = "exon";
        annotation.sjdbGTFtagExonParentTranscript = "transcript_id";
        annotation.sjdbGTFtagExonParentGene = "gene_id";
        annotation.sjdbGTFtagExonParentGeneName = {"gene_name"};
        annotation.sjdbGTFtagExonParentGeneType = {"gene_type", "gene_biotype"};
        annotation.sjdbOverhang = 100;
        annotation.limitSjdbInsertNsj = 1000000;
        std::uint64_t chromosome_bits = 18, prefix_bases = 14, runThreadN = 1;
        bool profile = false;
        auto last_stage = std::chrono::steady_clock::now();
        auto mark = [&](std::string_view name) {
            const auto now = std::chrono::steady_clock::now();
            if (profile)
                std::cerr << "GenomeGenerate stage " << name << ": "
                          << std::chrono::duration<double>(now-last_stage).count()
                          << " s\n";
            last_stage = now;
        };
        for (int i = first_option; i < argc;) {
            const std::string_view option(argv[i++]);
            auto value = [&]() -> std::string {
                if (i >= argc || std::string_view(argv[i]).starts_with("--"))
                    throw std::invalid_argument("missing value for " + std::string(option));
                return argv[i++];
            };
            if (option == "--runMode") {
                if (value() != "genomeGenerate")
                    throw std::invalid_argument("unsupported --runMode");
            } else if (option == "--profileGenomeGenerate") profile = true;
            else if (option == "--genomeFastaFiles") {
                const auto start = i;
                while (i < argc && !std::string_view(argv[i]).starts_with("--"))
                    build.fasta_files.emplace_back(argv[i++]);
                if (start == i) throw std::invalid_argument("--genomeFastaFiles needs a path");
            } else if (option == "--genomeDir") build.output_directory = value();
            else if (option == "--sjdbGTFfile") {
                const auto path = value();
                if (path == "-") annotation.sjdbGTFfile.reset();
                else annotation.sjdbGTFfile = path;
            }
            else if (option == "--sjdbOverhang") {
                annotation.sjdbOverhang = unsignedOption(value(), option, 1, 100000);
            } else if (option == "--genomeChrBinNbits")
                chromosome_bits = unsignedOption(value(), option, 8, 30);
            else if (option == "--runThreadN")
                runThreadN = unsignedOption(value(), "--runThreadN", 1, 256);
            else if (option == "--genomeSAindexNbases")
                prefix_bases = unsignedOption(value(), option, 1, 15);
            else if (option == "--limitSjdbInsertNsj")
                annotation.limitSjdbInsertNsj = unsignedOption(value(), option, 1, UINT64_MAX);
            else if (option == "--sjdbGTFchrPrefix") {
                const auto prefix = value();
                annotation.sjdbGTFchrPrefix = prefix == "-" ? "" : prefix;
            }
            else if (option == "--sjdbGTFfeatureExon") annotation.sjdbGTFfeatureExon = value();
            else if (option == "--sjdbGTFtagExonParentTranscript")
                annotation.sjdbGTFtagExonParentTranscript = value();
            else if (option == "--sjdbGTFtagExonParentGene")
                annotation.sjdbGTFtagExonParentGene = value();
            else if (option == "--sjdbGTFtagExonParentGeneName" ||
                     option == "--sjdbGTFtagExonParentGeneType") {
                auto& tags = option == "--sjdbGTFtagExonParentGeneName" ?
                    annotation.sjdbGTFtagExonParentGeneName : annotation.sjdbGTFtagExonParentGeneType;
                tags = {value()};
                while (i < argc && !std::string_view(argv[i]).starts_with("--"))
                    tags.emplace_back(argv[i++]);
            }
            else throw std::invalid_argument("unsupported genomeGenerate option: " +
                                             std::string(option));
        }
        if (build.fasta_files.empty() || build.output_directory.empty())
            throw std::invalid_argument("--genomeFastaFiles and --genomeDir are required");

        build.storage_options = BackendOptions{"nustar-reference-v1", 1,
            {OptionProvenance{"genomeChrBinNbits", {}, {std::to_string(chromosome_bits)}, "CLI"}}};
        GenomeReferenceService references;
        auto ready = take(references.buildBase(build));
        mark("FASTA/reference");
        if (annotation.sjdbGTFfile)
            ready = take(references.applyAnnotations(ready.lease,
                ReferenceAnnotationRequest{annotation, build.output_directory}));
        mark("GTF/SJDB");

        SearchIndexRequest search{"nustar-sa", build.output_directory / "index.sa",
            BackendOptions{"nustar-sa-v4", 4,
                {OptionProvenance{"genomeSAindexNbases", {},
                    {std::to_string(prefix_bases)}, "CLI"},
                 OptionProvenance{"runThreadN", {}, {std::to_string(runThreadN)}, "CLI"}}}};
        SuffixArrayIndexBackend backend;
        auto index = take(backend.build(search, ready.lease));
        mark("SAIS/prefix/write");

        // Verify both persisted artifacts before reporting a successful build.
        const auto saved_index_hash=index->identity().artifact_sha256;
        index.reset();
        ready.lease.reset();
        auto loaded_reference = take(references.loadBase(
            ReferenceRequest{build.output_directory, "nustar-reference-v1", {}}));
        auto loaded_index = take(backend.load(search, loaded_reference.lease));
        mark("reload/verify");
        const auto concrete = std::dynamic_pointer_cast<const GenomeReferenceData>(loaded_reference.lease);
        const auto suffix_index = std::dynamic_pointer_cast<const SuffixArrayIndex>(loaded_index);
        if (!concrete || !suffix_index || saved_index_hash !=
                loaded_index->identity().artifact_sha256)
            throw std::runtime_error("index reload verification failed");
        const auto manifest_path = build.output_directory / "nustarGenomeParameters.txt";
        std::ofstream manifest(manifest_path, std::ios::trunc);
        if (!manifest) throw std::runtime_error("cannot write build manifest");
        manifest << "builder\tnustar 0.2\n"
                 << "referenceFormat\tnustar-reference-v1\n"
                 << "indexFormat\tnustar-sa-v4\n"
                 << "saAlgorithm\t" << suffix_index->algorithmId() << '\n'
                 << "saPackedBitWidth\t"
                 << static_cast<unsigned>(suffix_index->suffixArray().bitWidth()) << '\n'
                 << "saPackedPayloadBytes\t"
                 << suffix_index->suffixArray().packedBytes().size() << '\n'
                 << "genomeChrBinNbits\t" << chromosome_bits << '\n'
                 << "genomeSAindexNbases\t" << prefix_bases << '\n'
                 << "sjdbOverhang\t" << (annotation.sjdbGTFfile ? annotation.sjdbOverhang : 0) << '\n'
                 << "referenceSemanticId\t" << concrete->semanticId() << '\n'
                 << "referenceSha256\t" << loaded_reference.files.front().sha256 << '\n'
                 << "indexSha256\t" << loaded_index->identity().artifact_sha256 << '\n';
        for (const auto& fasta : build.fasta_files)
            manifest << "genomeFastaFile\t" << fasta.string() << '\n';
        if (annotation.sjdbGTFfile)
            manifest << "sjdbGTFfile\t" << annotation.sjdbGTFfile->string() << '\n';
        manifest.close();
        if (!manifest) throw std::runtime_error("failed writing build manifest");
        std::cout << "GenomeGenerate complete: " << concrete->storage().metadata.chromosomes.size()
                  << " chromosomes, " << concrete->junctions().size() << " inserted junctions, "
                  << suffix_index->suffixArray().size() << " indexed suffixes\n"
                  << "Reference: " << (build.output_directory / "reference.bin") << '\n'
                  << "Suffix array index: " << search.artifact << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "GenomeGenerate: " << error.what() << '\n';
        return 2;
    }
}

int runUpgradeIndex(int argc, char** argv) {
    using namespace nustar::star2711b;
    if (argc!=4) {
        std::cerr << "Usage: nustar upgradeIndex SOURCE_DIR NEW_DIR\n";
        return 2;
    }
    try {
        const std::filesystem::path source(argv[2]), target(argv[3]);
        if (std::filesystem::exists(target))
            throw std::runtime_error("upgrade destination must not exist");
        auto reference=take(loadGenomeReference(source));
        SuffixArrayIndexBackend backend;
        auto index=take(backend.load({"nustar-sa",source/"index.sa",{}},reference));
        std::filesystem::create_directories(target);
        for (const auto& entry : std::filesystem::directory_iterator(source)) {
            const auto name=entry.path().filename();
            if (!entry.is_regular_file() || name=="index.sa" || name=="index.sa.sha256") continue;
            std::filesystem::copy_file(entry.path(),target/name);
        }
        auto saved=backend.save(*index,target/"index.sa");
        if (std::holds_alternative<Failure>(saved))
            throw std::runtime_error(std::get<Failure>(saved).message);
        // Keep the human-readable manifest consistent with the newly saved
        // artifact; source/reference provenance remains unchanged.
        std::ifstream checksum(target/"index.sa.sha256");
        std::string digest;
        if (!std::getline(checksum,digest)) throw std::runtime_error("missing upgraded checksum");
        const auto manifest_path=target/"nustarGenomeParameters.txt";
        std::ifstream previous(manifest_path);
        std::vector<std::string> lines;
        std::string line;
        bool format_seen=false, hash_seen=false;
        while (std::getline(previous,line)) {
            if (line.starts_with("indexFormat\t")) { line="indexFormat\tnustar-sa-v4"; format_seen=true; }
            if (line.starts_with("indexSha256\t")) { line="indexSha256\t"+digest; hash_seen=true; }
            lines.push_back(std::move(line));
        }
        previous.close();
        if (!format_seen) lines.push_back("indexFormat\tnustar-sa-v4");
        if (!hash_seen) lines.push_back("indexSha256\t"+digest);
        std::ofstream manifest(manifest_path,std::ios::trunc);
        for (const auto& row : lines) manifest << row << '\n';
        manifest.close();
        if (!manifest) throw std::runtime_error("cannot write upgraded manifest");
        std::cout << "Wrote packed multilevel prefix index v4 to " << target << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Index upgrade failed: " << e.what() << '\n';
        return 1;
    }
}
