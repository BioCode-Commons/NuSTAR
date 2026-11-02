#pragma once
#include "../output/sink.hpp"
#include "../reads/raw.hpp"
#include "config.hpp"
#include <functional>
namespace nustar::star2711b {
// The callback is used only after counting has finished, during sorted output.
using BamCorrection = std::function<std::pair<std::string, std::string>(std::uint64_t)>;
const IReadSourceBackend& bamSourceBackend();
Result<std::unique_ptr<IOutputSink>> openBamSink(const OutputFiles&, const BamConfig&,
                                                 const ReferenceMetadata&, BamCorrection = {});
} // namespace nustar::star2711b
