#include "genome.hpp"
#include "../common/binary.hpp"
#include "../common/file_identity.hpp"
#include "../common/sha256.hpp"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <fstream>
#include <limits>
#include <unordered_set>

namespace nustar::star2711b {
namespace {

std::atomic<std::uint64_t> next_generation{1};

Failure referenceError(ErrorKind kind, std::string message) {
    return Failure{kind, Stage::Reference, std::nullopt, 1, std::move(message)};
}

BaseCode decodeBase(char base) {
    switch (std::toupper(static_cast<unsigned char>(base))) {
        case 'A': return BaseCode::A;
        case 'C': return BaseCode::C;
        case 'G': return BaseCode::G;
        case 'T': case 'U': return BaseCode::T;
        default: return BaseCode::N;
    }
}

void hash64(Sha256& hash, std::uint64_t value) {
    char bytes[8];
    for (unsigned i = 0; i < 8; ++i) bytes[i] = char((value >> (8 * i)) & 255);
    hash.update(std::string_view(bytes, 8));
}

void writeJunction(std::ostream& out, const SjdbEntry& j) {
    using namespace binary;
    put64(out, j.id.value);
    put64(out, j.lookup_intron.begin.value); put64(out, j.lookup_intron.end.value);
    put64(out, j.reported_intron.begin.value); put64(out, j.reported_intron.end.value);
    out.put(char(j.strand)); out.put(char(j.motif));
    out.put(char(j.shift_left)); out.put(char(j.shift_right));
    put64(out, j.donor_sequence_begin.value); put64(out, j.acceptor_sequence_begin.value);
}

bool readJunction(std::istream& in, SjdbEntry& j) {
    using namespace binary;
    std::uint64_t values[5];
    for (auto& value : values) if (!get64(in, value)) return false;
    const int strand = in.get(), motif = in.get(), left = in.get(), right = in.get();
    if (strand < 0 || motif < 0 || left < 0 || right < 0 || strand > 2 || motif > 6) return false;
    std::uint64_t donor = 0, acceptor = 0;
    if (!get64(in, donor) || !get64(in, acceptor)) return false;
    j = SjdbEntry{JunctionId{values[0]},
        {GenomeOffset{values[1]}, GenomeOffset{values[2]}},
        {GenomeOffset{values[3]}, GenomeOffset{values[4]}},
        Strand(strand), Motif(motif), static_cast<std::uint8_t>(left),
        static_cast<std::uint8_t>(right), GenomeOffset{donor}, GenomeOffset{acceptor}};
    return true;
}

} // namespace

GenomeReferenceData::GenomeReferenceData(std::shared_ptr<const GenomeStorage> storage,
    std::vector<SjdbEntry> junctions, std::vector<BaseCode> synthetic,
    std::uint64_t sjdb_overhang, bool transient)
    : junctions_(std::move(junctions)), sjdb_overhang_(sjdb_overhang),
      generation_{next_generation.fetch_add(1)} {
    auto bytes=std::make_shared<std::vector<BaseCode>>();
    bytes->reserve(storage->forward.size()+synthetic.size()+40);
    bytes->insert(bytes->end(),20,BaseCode::GenomePadding);
    bytes->insert(bytes->end(),storage->forward.begin(),storage->forward.end());
    bytes->insert(bytes->end(),synthetic.begin(),synthetic.end());
    bytes->insert(bytes->end(),20,BaseCode::GenomePadding);
    auto owned=std::make_shared<GenomeStorage>(*storage);
    owned->sequence_storage=bytes;
    owned->forward=std::span<const BaseCode>(*bytes).subspan(20,storage->forward.size());
    storage_=std::move(owned);
    star_genome_=std::move(bytes);
    synthetic_=std::span<const BaseCode>(*star_genome_).subspan(
        20+storage_->forward.size(),synthetic.size());
    if(transient) semantic_id_="nustar-runtime-reference:"+std::to_string(generation_.value);
    else identify();
}

GenomeReferenceData::GenomeReferenceData(std::shared_ptr<const GenomeStorage> storage,
    std::vector<SjdbEntry> junctions, std::uint64_t overhang,
    std::shared_ptr<const std::vector<BaseCode>> genome, std::string semantic_id)
    : storage_(std::move(storage)), junctions_(std::move(junctions)),
      sjdb_overhang_(overhang), generation_{next_generation.fetch_add(1)},
      semantic_id_(std::move(semantic_id)),
      star_genome_(std::move(genome)) {
    synthetic_=std::span<const BaseCode>(*star_genome_).subspan(
        20+storage_->forward.size(),star_genome_->size()-40-storage_->forward.size());
}

void GenomeReferenceData::identify() {
    Sha256 hash;
    hash.update("nustar-reference-v1");
    hash64(hash, storage_->chromosome_bin_bits);
    for (const auto& chr : storage_->metadata.chromosomes) {
        hash.update(chr.name);
        hash64(hash, chr.start.value);
        hash64(hash, chr.length);
    }
    hash.update(std::as_bytes(std::span(storage_->forward)));
    hash64(hash, sjdb_overhang_);
    for (const auto& j : junctions_) {
        hash64(hash, j.lookup_intron.begin.value);
        hash64(hash, j.lookup_intron.end.value);
        hash64(hash, j.reported_intron.begin.value);
        hash64(hash, j.reported_intron.end.value);
        hash64(hash, std::uint64_t(j.strand));
        hash64(hash, std::uint64_t(j.motif));
        hash64(hash, j.shift_left);
        hash64(hash, j.shift_right);
    }
    hash.update(std::as_bytes(std::span(synthetic_)));
    semantic_id_ = hash.finish();
}

ReferenceView GenomeReferenceData::view() const {
    return ReferenceView{generation_, semantic_id_, &storage_->metadata, this,
                         JunctionDatabaseView{junctions_, sjdb_overhang_}};
}

std::span<const BaseCode> GenomeReferenceData::starGenome() const {
    return std::span<const BaseCode>(*star_genome_).subspan(20,star_genome_->size()-40);
}

Result<SequenceWindow> GenomeReferenceData::pin(GenomeRange range) const {
    return SequenceWindow{range,
        std::span(storage_->forward).subspan(range.begin.value, range.end.value - range.begin.value),
        std::static_pointer_cast<const ISequenceLease>(shared_from_this())};
}

Result<std::shared_ptr<const GenomeStorage>> readGenomeFasta(
    std::span<const std::filesystem::path> files, std::uint32_t chromosome_bin_bits) {
    if (files.empty() || chromosome_bin_bits < 8 || chromosome_bin_bits > 30) {
        return referenceError(ErrorKind::Config, "FASTA files or chromosome bin bits are invalid");
    }
    auto storage = std::make_shared<GenomeStorage>();
    auto sequence=std::make_shared<std::vector<BaseCode>>();
    auto& forward=*sequence;
    storage->chromosome_bin_bits = chromosome_bin_bits;
    storage->metadata.coordinate_schema = "star2711b-padded-bin-" + std::to_string(chromosome_bin_bits);
    const std::uint64_t bin = std::uint64_t{1} << chromosome_bin_bits;
    std::unordered_set<std::string> seen;
    bool have_chr = false;

    auto finish_chromosome = [&]() -> Status {
        if (!have_chr) return std::monostate{};
        auto& chr = storage->metadata.chromosomes.back();
        chr.length = forward.size() - chr.start.value;
        if (chr.length == 0) return referenceError(ErrorKind::Input, "empty FASTA chromosome: " + chr.name);
        const auto n = forward.size();
        // STAR genomeScanFastaFiles leaves an extra bin at a boundary.
        const auto padded = ((n + 1) / bin + 1) * bin;
        forward.resize(padded, BaseCode::GenomePadding);
        return std::monostate{};
    };

    for (const auto& file : files) {
        std::ifstream input(file);
        if (!input) return referenceError(ErrorKind::Io, "cannot open FASTA: " + file.string());
        std::string line;
        while (std::getline(input, line)) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                if (auto status = finish_chromosome(); std::holds_alternative<Failure>(status))
                    return std::get<Failure>(std::move(status));
                const auto end = line.find_first_of(" \t\r", 1);
                const auto name = line.substr(1, end == std::string::npos ? end : end - 1);
                if (name.empty() || !seen.insert(name).second) {
                    return referenceError(ErrorKind::Input, "empty or duplicate FASTA chromosome: " + name);
                }
                storage->metadata.chromosomes.push_back(Chromosome{
                    ChromosomeId{storage->metadata.chromosomes.size()}, name,
                    GenomeOffset{forward.size()}, 0});
                have_chr = true;
            } else {
                if (!have_chr) return referenceError(ErrorKind::Input, "FASTA sequence precedes a header");
                for (unsigned char base : line) {
                    if (!std::isspace(base)) forward.push_back(decodeBase(base));
                }
            }
        }
        if (!input.eof()) return referenceError(ErrorKind::Io, "error reading FASTA: " + file.string());
    }
    if (!have_chr) return referenceError(ErrorKind::Input, "FASTA contains no chromosomes");
    if (auto status = finish_chromosome(); std::holds_alternative<Failure>(status))
        return std::get<Failure>(std::move(status));
    storage->metadata.real_genome_end = GenomeOffset{forward.size()};
    storage->sequence_storage=sequence;
    storage->forward=*sequence;
    return std::shared_ptr<const GenomeStorage>(std::move(storage));
}

Result<FileIdentity> saveGenomeReference(const GenomeReferenceData& reference,
    const std::filesystem::path& directory) {
    std::error_code error;
    std::filesystem::create_directories(directory, error);
    if (error) return referenceError(ErrorKind::Io, "cannot create genome directory: " + error.message());
    const auto path = directory / "reference.bin";
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) return referenceError(ErrorKind::Io, "cannot create reference: " + path.string());
    out.write("KREF", 4);
    binary::put32(out, 1);
    binary::put32(out, reference.storage().chromosome_bin_bits);
    binary::put32(out, static_cast<std::uint32_t>(reference.storage().metadata.chromosomes.size()));
    binary::put64(out, reference.storage().forward.size());
    binary::put64(out, reference.junctions().size());
    binary::put64(out, reference.synthetic().size());
    binary::put64(out, reference.sjdbOverhang());
    for (const auto& chr : reference.storage().metadata.chromosomes) {
        binary::putString(out, chr.name);
        binary::put64(out, chr.start.value);
        binary::put64(out, chr.length);
    }
    out.write(reinterpret_cast<const char*>(reference.storage().forward.data()),
              static_cast<std::streamsize>(reference.storage().forward.size()));
    for (const auto& j : reference.junctions()) writeJunction(out, j);
    out.write(reinterpret_cast<const char*>(reference.synthetic().data()),
              static_cast<std::streamsize>(reference.synthetic().size()));
    binary::putString(out, std::string(reference.semanticId()));
    out.close();
    if (!out) return referenceError(ErrorKind::Io, "failed writing reference: " + path.string());
    std::ofstream names(directory / "chrName.txt");
    std::ofstream lengths(directory / "chrLength.txt");
    std::ofstream starts(directory / "chrStart.txt");
    std::ofstream name_lengths(directory / "chrNameLength.txt");
    if (!names || !lengths || !starts || !name_lengths)
        return referenceError(ErrorKind::Io, "cannot create chromosome metadata tables");
    for (const auto& chr : reference.storage().metadata.chromosomes) {
        names << chr.name << '\n';
        lengths << chr.length << '\n';
        starts << chr.start.value << '\n';
        name_lengths << chr.name << '\t' << chr.length << '\n';
    }
    starts << reference.storage().forward.size() << '\n';
    names.close(); lengths.close(); starts.close(); name_lengths.close();
    if (!names || !lengths || !starts || !name_lengths)
        return referenceError(ErrorKind::Io, "failed writing chromosome metadata tables");
    return identifyFile(path);
}

Result<std::shared_ptr<const GenomeReferenceData>> loadGenomeReference(
    const std::filesystem::path& directory) {
    const auto path = directory / "reference.bin";
    std::ifstream in(path, std::ios::binary);
    if (!in) return referenceError(ErrorKind::Io, "cannot open reference: " + path.string());
    char magic[4]{};
    in.read(magic, 4);
    std::uint32_t version = 0, bin_bits = 0, chr_count = 0;
    std::uint64_t forward_size = 0, sj_count = 0, synthetic_size = 0, overhang = 0;
    if (std::string_view(magic, 4) != "KREF" ||
        !binary::get32(in, version) || version != 1 ||
        !binary::get32(in, bin_bits) || bin_bits < 8 || bin_bits > 30 ||
        !binary::get32(in, chr_count) || chr_count == 0 || chr_count > 100000 ||
        !binary::get64(in, forward_size) || !binary::get64(in, sj_count) ||
        !binary::get64(in, synthetic_size) || !binary::get64(in, overhang)) {
        return referenceError(ErrorKind::Index, "invalid reference header: " + path.string());
    }
    std::error_code ec;
    const auto file_size = std::filesystem::file_size(path, ec);
    if (ec || forward_size > file_size || synthetic_size > file_size || sj_count > file_size / 60) {
        return referenceError(ErrorKind::Index, "invalid reference lengths: " + path.string());
    }
    auto storage = std::make_shared<GenomeStorage>();
    storage->chromosome_bin_bits = bin_bits;
    storage->metadata.coordinate_schema = "star2711b-padded-bin-" + std::to_string(bin_bits);
    for (std::uint32_t i = 0; i < chr_count; ++i) {
        std::string name;
        std::uint64_t start = 0, length = 0;
        if (!binary::getString(in, name) || !binary::get64(in, start) ||
            !binary::get64(in, length) || start > forward_size || length > forward_size - start)
            return referenceError(ErrorKind::Index, "invalid reference chromosome table");
        storage->metadata.chromosomes.push_back(Chromosome{ChromosomeId{i},
            std::move(name), GenomeOffset{start}, length});
    }
    if (forward_size>std::numeric_limits<std::size_t>::max()-40 ||
        synthetic_size>std::numeric_limits<std::size_t>::max()-40-forward_size)
        return referenceError(ErrorKind::Capacity,"reference buffer length overflow");
    auto sequence=std::make_shared<std::vector<BaseCode>>(
        forward_size+synthetic_size+40,BaseCode::GenomePadding);
    storage->sequence_storage=sequence;
    storage->forward=std::span<const BaseCode>(*sequence).subspan(20,forward_size);
    in.read(reinterpret_cast<char*>(sequence->data()+20),
            static_cast<std::streamsize>(forward_size));
    std::vector<SjdbEntry> junctions(static_cast<std::size_t>(sj_count));
    for (auto& j : junctions) if (!readJunction(in, j))
        return referenceError(ErrorKind::Index, "invalid junction table");
    in.read(reinterpret_cast<char*>(sequence->data()+20+forward_size),
            static_cast<std::streamsize>(synthetic_size));
    std::string stored_semantic;
    if (!in || !binary::getString(in, stored_semantic, 128))
        return referenceError(ErrorKind::Index, "truncated reference file");
    storage->metadata.real_genome_end = GenomeOffset{forward_size};
    auto reference=std::shared_ptr<const GenomeReferenceData>(new GenomeReferenceData(
        std::move(storage),std::move(junctions),overhang,std::move(sequence),
        std::move(stored_semantic)));
    return reference;
}

} // namespace nustar::star2711b
