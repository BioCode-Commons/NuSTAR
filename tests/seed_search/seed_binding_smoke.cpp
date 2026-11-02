#include "../../runtime/reference_service.hpp"
#include "../../seed_search/star_search.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <filesystem>
#include <exception>
#include <iostream>
#include <random>
#include <sstream>
#include <thread>
#include <stdexcept>

using namespace nustar::star2711b;

namespace {
template<class T> T checked(Result<T> result) {
    if (std::holds_alternative<Failure>(result))
        throw std::runtime_error(std::get<Failure>(result).message);
    return std::get<T>(std::move(result));
}
void checked(Status result) {
    if (std::holds_alternative<Failure>(result))
        throw std::runtime_error(std::get<Failure>(result).message);
}
std::string snapshot(const StarSeedSearchResult& result) {
    std::ostringstream out;
    const auto& c = result.counters;
    out << result.map_marker << ':' << c.all_alignments << ':' << c.unique_alignments
        << ':' << c.multiple_alignments << ':' << c.min_multimaps << ':' << c.min_multimap_length
        << ':' << c.max_unique_length << ':' << c.no_good_piece_length;
    for (const auto& split : result.splits)
        out << '|' << split.read.begin.value << ',' << split.read.end.value << ',' << int(split.mate);
    for (const auto& p : result.pieces)
        out << ';' << p.read_start << ',' << p.length << ',' << p.direction << ',' << p.multiplicity
            << ',' << p.sa_begin << ',' << p.sa_end << ',' << p.fragment;
    return out.str();
}
}

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    const auto directory = std::filesystem::temp_directory_path() /
        ("nustar-seed-smoke-" + std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count()));
    struct Cleanup {
        std::filesystem::path path;
        ~Cleanup() { std::error_code error; std::filesystem::remove_all(path, error); }
    } cleanup{directory};
    try {
        GenomeReferenceService references;
        auto ready = checked(references.buildBase(ReferenceBuildRequest{
            {std::filesystem::path(argv[1])}, directory,
            BackendOptions{"nustar-reference-v1", 1,
                {OptionProvenance{"genomeChrBinNbits", {}, {"8"}, "test"}}}}));
        AnnotationConfig annotation{};
        annotation.sjdbGTFfile = std::filesystem::path(argv[2]);
        annotation.sjdbGTFfeatureExon = "exon";
        annotation.sjdbGTFtagExonParentTranscript = "transcript_id";
        annotation.sjdbGTFtagExonParentGene = "gene_id";
        annotation.sjdbOverhang = 16;
        annotation.limitSjdbInsertNsj = 100;
        ready = checked(references.applyAnnotations(
            ready.lease, ReferenceAnnotationRequest{annotation, directory}));
        auto reference = std::dynamic_pointer_cast<const GenomeReferenceData>(ready.lease);
        SuffixArrayIndexBackend backend;
        SearchIndexRequest request{"nustar-sa", directory/"index.sa",
            BackendOptions{"nustar-sa-v3", 3,
                {OptionProvenance{"genomeSAindexNbases", {}, {"4"}, "test"}}}};
        auto index = std::dynamic_pointer_cast<const SuffixArrayIndex>(
            checked(backend.build(request, reference)));
        if (!index || reference->junctions().empty())
            throw std::runtime_error("test index is incomplete");
        MappingRead read{};
        read.fragment = FragmentId{0, 17};
        read.encoded[0].assign(reference->synthetic().begin()+8,
                               reference->synthetic().begin()+24);
        read.encoded[1].reserve(read.encoded[0].size());
        for (auto base : read.encoded[0]) {
            const auto code = static_cast<std::uint8_t>(base);
            read.encoded[1].push_back(code < 4 ?
                static_cast<BaseCode>(code^3) : base);
        }
        read.encoded[2].assign(read.encoded[1].rbegin(), read.encoded[1].rend());
        read.mates.push_back(PreparedMate{16, 16, {},
            ReadRange{ReadOffset{0}, ReadOffset{16}}, false});
        read.budget = ReadBudget{0, 16, 16, 16};

        const SeedConfig config{50, 0, 1.0, 10000, 1000, 12, 5, 10, 0};
        StarSeedSearch native(index, reference, config);
        const auto native_result = checked(native.search(read));
        SeedDeliveryContract delivery{SeedContract::Star2711bEquivalent,
            HitOrder::Star2711b, true, true, true, 1,
            std::string(star_seed_search_domain)};
        auto plan = std::make_shared<SeedSearchPlan>(SeedSearchPlan{
            delivery, std::nullopt, config, BackendOptions{}});
        StarSeedAlgorithmFactory factory;
        auto finder = checked(factory.bind(index, reference, plan));
        auto workspace = checked(finder->createWorkspace({}, WorkspaceBudget{}));
        SeedBatch batch{};
        checked(finder->find(read, *workspace, batch));
        if (batch.seeds.size() != native_result.pieces.size() ||
            batch.splits.size() != native_result.splits.size() ||
            batch.markers.legacy_map_marker != native_result.map_marker ||
            !batch.occurrences)
            throw std::runtime_error("public seed output differs from native search");
        bool crossing_junction = false;
        for (std::size_t i = 0; i < batch.seeds.size(); ++i) {
            const auto& seed = batch.seeds[i];
            const auto& piece = native_result.pieces[i];
            if (seed.read.begin.value != piece.read_start ||
                seed.read.end.value != piece.read_start+piece.length ||
                seed.multiplicity.value != piece.multiplicity ||
                seed.multiplicity.precision != CountPrecision::Exact)
                throw std::runtime_error("public seed fields differ from STAR PC");
            auto cursor = checked(batch.occurrences->open(seed.occurrences));
            std::uint64_t count = 0;
            SeedHitPage page{};
            do {
                checked(cursor->next(3, page));
                for (const auto& hit : page.hits) {
                    if (hit.visit_ordinal != count++)
                        throw std::runtime_error("seed hit order changed");
                    if (hit.blocks.size() == 2 &&
                        hit.blocks[0].junction && hit.blocks[1].junction)
                        crossing_junction = true;
                }
            } while (!page.end);
            // STAR counts every matching SA suffix in multiplicity, but its
            // junction projection discards synthetic hits that do not cross
            // the donor/acceptor boundary.
            if (count > piece.multiplicity)
                throw std::runtime_error("seed cursor invented occurrences");
            std::uint64_t visited = 0;
            checked(batch.occurrences->visit(seed.occurrences,
                [&](const SeedHitView& hit) {
                    if (hit.visit_ordinal != visited++)
                        throw std::runtime_error("seed visitor order changed");
                    return true;
                }));
            if (visited != count)
                throw std::runtime_error("seed cursor and visitor disagree");
            visited=0;
            checked(batch.occurrences->visitBatch(seed.occurrences,
                [&](std::span<const SeedHitView> hits) {
                    for (const auto& hit : hits)
                        if (hit.visit_ordinal!=visited++)
                            throw std::runtime_error("batched seed order changed");
                    return true;
                }));
            if (visited!=count) throw std::runtime_error("batched seed count differs");
            unsigned calls=0;
            checked(batch.occurrences->visitBatch(seed.occurrences,
                [&](std::span<const SeedHitView> hits) {
                    if (hits.empty() || calls++) throw std::runtime_error("batch early stop ignored");
                    return false;
                }));
            if (calls!=unsigned(count!=0)) throw std::runtime_error("batch early stop lost hits");
        }
        // The light anchor pass must match full projection, including donor/
        // acceptor order and stopping between the two sides of a junction.
        for (const auto& seed : batch.seeds) {
            std::vector<std::pair<std::uint64_t, unsigned>> full, anchors;
            checked(batch.occurrences->visit(seed.occurrences, [&](const SeedHitView& hit) {
                for (const auto& b : hit.blocks)
                    full.emplace_back(b.genome_start.value, unsigned(hit.orientation));
                return true;
            }));
            std::vector<std::string> complete_blocks, direct_blocks;
            const auto block_text = [](const SeedHitBlock& block, Orientation orientation) {
                std::ostringstream out;
                out << block.read_start.value << ',' << block.genome_start.value << ','
                    << block.length << ',' << unsigned(block.mate) << ','
                    << (block.junction ? block.junction->value : UINT64_MAX) << ','
                    << unsigned(orientation);
                return out.str();
            };
            checked(batch.occurrences->visit(seed.occurrences, [&](const SeedHitView& hit) {
                for (const auto& block : hit.blocks) complete_blocks.push_back(block_text(block, hit.orientation));
                return true;
            }));
            checked(batch.occurrences->visitBlocks(seed.occurrences,
                [&](const SeedHitBlock& block, Orientation orientation) {
                    direct_blocks.push_back(block_text(block, orientation)); return true;
                }));
            if (direct_blocks != complete_blocks) throw std::runtime_error("direct block projection differs");
            unsigned direct_calls = 0;
            checked(batch.occurrences->visitBlocks(seed.occurrences,
                [&](const SeedHitBlock&, Orientation) { ++direct_calls; return false; }));
            if (direct_calls != unsigned(!complete_blocks.empty()))
                throw std::runtime_error("direct block early stop differs");
            checked(batch.occurrences->visitAnchors(seed.occurrences,
                [&](GenomeOffset start, Orientation direction) {
                    anchors.emplace_back(start.value, unsigned(direction)); return true;
                }));
            if (anchors != full) throw std::runtime_error("anchor projection changed order/coordinates");
            for (std::size_t stop = 1; stop <= full.size(); ++stop) {
                anchors.clear();
                checked(batch.occurrences->visitAnchors(seed.occurrences,
                    [&](GenomeOffset start, Orientation direction) {
                        anchors.emplace_back(start.value, unsigned(direction));
                        return anchors.size() < stop;
                    }));
                if (anchors != decltype(full)(full.begin(), full.begin()+stop))
                    throw std::runtime_error("anchor early stop changed");
            }
        }
        // A retained source or cursor must not be recycled by the next find().
        auto retained = batch;
        const auto retained_id = retained.occurrences->identity().source;
        auto cursor = checked(retained.occurrences->open(retained.seeds.front().occurrences));
        for (unsigned iteration = 0; iteration < 3; ++iteration) {
            checked(finder->find(read, *workspace, batch));
            if (batch.occurrences->identity().source.value == retained_id.value ||
                retained.occurrences->identity().source.value != retained_id.value)
                throw std::runtime_error("retained seed source was mutated");
        }
        SeedHitPage retained_page{};
        checked(cursor->next(10000, retained_page));
        if (retained_page.hits.empty()) throw std::runtime_error("retained cursor lost hits");
        if (!crossing_junction)
            throw std::runtime_error("junction-spanning occurrence was not projected");
        // Recycled storage must not revive an old weak lease (or mutate an
        // observer which locks it concurrently). Releasing the last published
        // owner ends that control block, even while its storage is cached.
        std::weak_ptr<const ISeedOccurrenceSource> old_weak = batch.occurrences;
        batch.occurrences.reset();
        if (!old_weak.expired()) throw std::runtime_error("unpublished source kept a weak lease alive");
        checked(finder->find(read, *workspace, batch));
        if (!old_weak.expired()) throw std::runtime_error("source reuse revived a weak lease");
        // Occurrence ownership must survive workspace destruction. Concurrent
        // workers must never generate aliasing source handles, even across ID
        // allocation blocks and repeated searches for the same fragment.
        workspace.reset();
        checked(batch.occurrences->visit(batch.seeds.front().occurrences,
            [](const SeedHitView&) { return true; }));
        std::array<std::vector<std::uint64_t>, 8> ids;
        std::array<std::exception_ptr, 8> errors{};
        std::vector<std::jthread> workers;
        for (std::size_t worker = 0; worker < ids.size(); ++worker)
            workers.emplace_back([&, worker] {
                try {
                    auto scratch = checked(finder->createWorkspace({}, {}));
                    for (unsigned n = 0; n < 4200; ++n) {
                        SeedBatch current;
                        checked(finder->find(read, *scratch, current));
                        ids[worker].push_back(current.occurrences->identity().source.value);
                    }
                } catch (...) { errors[worker] = std::current_exception(); }
            });
        for (auto& worker : workers) worker.join();
        std::vector<std::uint64_t> all_ids;
        for (std::size_t worker = 0; worker < ids.size(); ++worker) {
            if (errors[worker]) std::rethrow_exception(errors[worker]);
            all_ids.insert(all_ids.end(), ids[worker].begin(), ids[worker].end());
        }
        std::sort(all_ids.begin(), all_ids.end());
        if (std::adjacent_find(all_ids.begin(), all_ids.end()) != all_ids.end())
            throw std::runtime_error("concurrent seed source identities collided");
        // A sparse 14-base directory must return exactly the same ordered PC
        // and counters as the independently bounded dense 4-base search.
        request.artifact = directory/"index14.sa";
        request.options.values[0].resolved_tokens = {"14"};
        auto wide_index = std::dynamic_pointer_cast<const SuffixArrayIndex>(
            checked(backend.build(request, reference)));
        StarSeedSearch sparse(wide_index, reference, config);
        std::mt19937 random(991);
        const auto& forward = reference->storage().forward;
        for (unsigned trial = 0; trial < 500; ++trial) {
            const auto length = 5+trial%60;
            const auto start = random()%(forward.size()-length);
            read.encoded[0].assign(forward.begin()+start, forward.begin()+start+length);
            if (trial%3 == 0)
                for (auto& b : read.encoded[0]) b = BaseCode(random()%5);
            if (trial%13 == 0)
                std::fill(read.encoded[0].begin(), read.encoded[0].end(), BaseCode((trial/13)%4));
            read.encoded[1].clear();
            for (auto base : read.encoded[0]) {
                const auto code = static_cast<unsigned>(base);
                read.encoded[1].push_back(code < 4 ? BaseCode(code^3) : base);
            }
            read.encoded[2].clear();
            for (auto it = read.encoded[0].rbegin(); it != read.encoded[0].rend(); ++it) {
                const auto code = static_cast<unsigned>(*it);
                read.encoded[2].push_back(code < 4 ? BaseCode(code^3) : *it);
            }
            read.mates[0] = PreparedMate{length, length, {}, {{0}, {length}}, false};
            read.budget = ReadBudget{0, length, length, length};
            if (snapshot(checked(native.search(read))) != snapshot(checked(sparse.search(read))))
                throw std::runtime_error("dense/sparse prefix lookup changed seeds or SA rank intervals");
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
