#include "../../reads/raw.hpp"

#include <cassert>
#include <utility>
#include <variant>

using namespace nustar::star2711b;

int main() {
    ReadBatch se{BatchId{7}, {}};
    se.fragments.push_back(RawFragment{
        FragmentId{0, 9}, Layout::SE, FastxFormat::Fasta,
        {RawMate{"read", "", std::nullopt, ""}, std::nullopt}, 0});
    auto single = takeSingleBatch(std::move(se));
    assert(std::holds_alternative<SingleReadBatch>(single));
    assert(std::get<SingleReadBatch>(single).reads[0].mate.sequence.empty());

    ReadBatch pe{BatchId{8}, {}};
    pe.fragments.push_back(RawFragment{
        FragmentId{1, 10}, Layout::PE, FastxFormat::Fastq,
        {RawMate{"first", "A", std::string("!"), ""}, std::nullopt}, 0});
    auto incomplete = takePairedBatch(std::move(pe));
    assert(std::holds_alternative<Failure>(incomplete));
    assert(std::get<Failure>(incomplete).stage == Stage::Input);
    assert(std::get<Failure>(incomplete).fragment->ordinal == 10);
}
