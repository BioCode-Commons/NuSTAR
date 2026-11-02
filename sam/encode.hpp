#pragma once

#include "model.hpp"
#include "config.hpp"
#include "../reference/view.hpp"

namespace nustar::star2711b {

struct EncodedFragment {
    FragmentId fragment;
    std::vector<std::byte> sam;
    std::array<std::vector<std::byte>, 2> unmapped_fastx;
};
Result<EncodedFragment> encodeSam(const FragmentId&, const ProjectedOutput&);
Result<std::vector<std::byte>> encodeSamHeader(const ReferenceMetadata&,
    const SamConfig&, const BaselineIdentity&);

} // namespace nustar::star2711b
