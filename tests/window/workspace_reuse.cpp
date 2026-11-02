#include "../../window/star_window.hpp"

#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>

using namespace nustar::star2711b;

namespace {
template<class T> T checked(Result<T> result) {
    if (auto* error = std::get_if<Failure>(&result)) throw std::runtime_error(error->message);
    return std::get<T>(std::move(result));
}
void checked(Status result) {
    if (auto* error = std::get_if<Failure>(&result)) throw std::runtime_error(error->message);
}
struct Reference final : IReferenceData {
    ReferenceMetadata metadata{{{ChromosomeId{0}, "a", GenomeOffset{0}, 512},
        {ChromosomeId{1}, "b", GenomeOffset{512}, 512}}, GenomeOffset{1024}, "test"};
    ReferenceView view() const override { return {{1}, "test", &metadata, nullptr, {}}; }
};
struct Source final : ISeedOccurrenceSource {
    OccurrenceSourceIdentity id{{1}, {0, 0}, {1}, "test",
        {SeedContract::Star2711bEquivalent, HitOrder::Star2711b, true, true, true, 1, "test"}};
    std::vector<SeedHit> hits;
    bool fail = false;
    const OccurrenceSourceIdentity& identity() const noexcept override { return id; }
    Result<std::unique_ptr<ISeedHitCursor>> open(OccurrenceHandle) const override {
        throw std::logic_error("visitor expected");
    }
    Status visit(OccurrenceHandle, const std::function<bool(const SeedHitView&)>& fn) const override {
        for (const auto& hit : hits) {
            if (!fn({hit.chromosome, hit.orientation, hit.visit_ordinal, hit.blocks})) break;
            if (fail) return Failure{ErrorKind::Input, Stage::Window, {}, 1, "injected failure"};
        }
        return std::monostate{};
    }
};
std::string snapshot(const WindowBatch& batch) {
    std::ostringstream out;
    out << batch.generation.value << ':' << batch.markers.normal_early_stop
        << ':' << batch.markers.legacy_map_marker;
    for (const auto& w : batch.windows) {
        out << '|' << w.id.value << ',' << w.chromosome.value << ',' << int(w.orientation)
            << ',' << w.first_bin.value << ',' << w.last_bin_inclusive.value << ',' << w.alive
            << ',' << w.minimum_recordable_length << ',' << w.last_anchor.value_or(UINT64_MAX);
        for (const auto& s : w.seeds)
            out << ';' << s.read_start.value << ',' << s.genome_start.value << ',' << s.length
                << ',' << s.repeat_count << ',' << s.anchor << ',' << int(s.mate)
                << ',' << (s.sjdb ? s.sjdb->value : UINT64_MAX);
    }
    return out.str();
}
}
int main() {
    try {
        auto reference = std::make_shared<Reference>();
        auto builder = checked(makeStarWindowBuilder(reference));
        auto source = std::make_shared<Source>();
        MappingRead read{};
        read.budget.legacy_Lread = 100;
        SeedBatch seeds{};
        seeds.generation = {1}; seeds.delivery = source->id.delivery; seeds.occurrences = source;
        seeds.seeds.push_back({{{0}, {20}}, MateId::First, {1, CountPrecision::Exact}, {{1}, {0}}});
        std::mt19937 random(31);
        // Include capacity rollback, seed capacity early stop, both strands,
        // chromosome boundaries, merges and overlapping flanks.
        for (const auto capacity : {2U, 100U}) {
            WindowConfig config{50, 1, 0, 512, 9, 4, 4, 10, capacity, 0, 0};
            auto reused = checked(builder->createWorkspace(config, reference->view(), {}, {}));
            WindowBatch actual;
            for (unsigned iteration = 0; iteration < 400; ++iteration) {
                source->hits.clear();
                const auto count = iteration%13;
                for (unsigned i = 0; i < count; ++i) {
                    const auto bin = random()%512;
                    source->hits.push_back({ChromosomeId{bin < 256 ? 0U : 1U},
                        Orientation(i%2), i, {{{i*3}, {bin*2}, 20, MateId::First, {}}}});
                }
                // Also exercise a failed build after mutating the map.
                source->fail = iteration%11 == 0;
                WindowBatch unused;
                builder->build(read, seeds, reference->view(), config, *reused, unused);
                source->fail = false;
                if (iteration%3 == 0) { reused->reset(); reused->reset(); }
                seeds.markers.normal_early_stop = iteration%7 == 0;
                auto fresh = checked(builder->createWorkspace(config, reference->view(), {}, {}));
                WindowBatch expected;
                checked(builder->build(read, seeds, reference->view(), config, *reused, actual));
                checked(builder->build(read, seeds, reference->view(), config, *fresh, expected));
                if (snapshot(actual) != snapshot(expected))
                    throw std::runtime_error("window workspace retained previous read state");
            }
        }
        // STAR replacement ordering: move left/right/in place, and insert
        // after existing equal starts. The longer seed must actually overlap
        // under STAR's endpoint test (containing an old seed is not sufficient).
        seeds.markers = {};
        WindowConfig config{50, 9, 0, 2, 9, 4, 4, 100, 100, 0, 0};
        auto workspace = checked(builder->createWorkspace(config, reference->view(), {}, {}));
        auto replacement = [&](std::initializer_list<std::array<std::uint64_t, 3>> input,
                               std::initializer_list<std::array<std::uint64_t, 3>> expected) {
            source->hits.clear();
            for (const auto& value : input)
                source->hits.push_back({{0}, Orientation::Forward, source->hits.size(),
                    {{{value[0]}, {value[1]}, value[2], MateId::First, {}}}});
            WindowBatch result;
            checked(builder->build(read, seeds, reference->view(), config, *workspace, result));
            std::vector<std::array<std::uint64_t, 3>> actual;
            for (const auto& window : result.windows)
                for (const auto& seed : window.seeds)
                    actual.push_back({seed.read_start.value, seed.genome_start.value, seed.length});
            if (actual != std::vector<std::array<std::uint64_t, 3>>(expected))
                throw std::runtime_error("STAR overlapping seed replacement order changed");
        };
        replacement({{30,130,10},{20,220,5},{15,115,20}}, {{15,115,20},{20,220,5}});
        replacement({{10,110,30},{15,215,5},{20,120,40}}, {{15,215,5},{20,120,40}});
        replacement({{20,120,10},{30,230,5},{20,120,20}}, {{20,120,20},{30,230,5}});
        replacement({{25,125,10},{20,220,5},{20,320,5},{20,120,12}},
                    {{20,220,5},{20,320,5},{20,120,12}});
        return 0;
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
