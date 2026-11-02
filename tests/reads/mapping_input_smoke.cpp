#include "../../reads/model.hpp"

#include <cassert>
#include <utility>
#include <variant>

using namespace nustar::star2711b;

int main() {
    MappingRead single{};
    single.fragment = FragmentId{0, 1};
    single.mates.push_back(PreparedMate{2, 2, {}, {ReadOffset{0}, ReadOffset{2}}, false});
    for (auto& encoding : single.encoded) {
        encoding = {BaseCode::A, BaseCode::C};
    }
    single.budget.legacy_Lread = 5; // PE merged read may inherit a larger budget.
    auto accepted = SingleMappingRead::fromPrepared(std::move(single));
    assert(accepted.view().mates.size() == 1);
    assert(accepted.view().budget.legacy_Lread == 5);
    assert(accepted.view().encoded[0].size() == 2);

    MappingRead pair{};
    pair.fragment = FragmentId{0, 2};
    pair.mates.push_back(PreparedMate{1, 1, {}, {ReadOffset{0}, ReadOffset{1}}, false});
    pair.mates.push_back(PreparedMate{1, 1, {}, {ReadOffset{2}, ReadOffset{3}}, true});
    for (auto& encoding : pair.encoded) {
        encoding = {BaseCode::A, BaseCode::MateSpacer, BaseCode::T};
    }
    auto paired = PairedMappingRead::fromPrepared(std::move(pair));
    assert(paired.view().mates.size() == 2);
    assert(paired.view().encoded[0][1] == BaseCode::MateSpacer);
    assert(paired.view().mates[1].combined_forward_interval.begin.value == 2);
}
