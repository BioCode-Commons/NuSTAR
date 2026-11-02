#include "raw.hpp"

#include <utility>

namespace nustar::star2711b {
namespace {

Failure invalidFragment(FragmentId id, std::string message) {
    return Failure{ErrorKind::Input, Stage::Input, id, 1, std::move(message)};
}

} // namespace

Result<SingleReadBatch> takeSingleBatch(ReadBatch&& input) {
    SingleReadBatch output{input.id, {}};
    output.reads.reserve(input.fragments.size());
    for (auto& fragment : input.fragments) {
        if (fragment.layout != Layout::SE || !fragment.mates[0] || fragment.mates[1]) {
            return invalidFragment(fragment.id, "SE batch requires mate 1 and no mate 2");
        }
        output.reads.push_back(RawSingleRead{
            fragment.id, fragment.format, std::move(*fragment.mates[0]),
            fragment.legacy_read_filter});
    }
    return output;
}

Result<PairedReadBatch> takePairedBatch(ReadBatch&& input) {
    PairedReadBatch output{input.id, {}};
    output.pairs.reserve(input.fragments.size());
    for (auto& fragment : input.fragments) {
        if (fragment.layout != Layout::PE || !fragment.mates[0] || !fragment.mates[1]) {
            return invalidFragment(fragment.id, "PE batch requires both mates");
        }
        output.pairs.push_back(RawReadPair{
            fragment.id, fragment.format,
            {std::move(*fragment.mates[0]), std::move(*fragment.mates[1])},
            fragment.legacy_read_filter});
    }
    return output;
}

} // namespace nustar::star2711b
