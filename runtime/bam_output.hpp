#pragma once
#include "../bam/io.hpp"
#include "../solo/session.hpp"
namespace nustar::star2711b {
inline Result<std::unique_ptr<IOutputSink>>
openAlignmentSink(const OutputFiles& files, const BamConfig& config, const ReferenceMetadata& ref,
                  const std::shared_ptr<solo::Session>& solo) {
    if (!config.enabled())
        return bufferedSamSinkBackend().open(files);
    BamCorrection correction;
    if (solo && solo->hasCorrections())
        correction = [solo](std::uint64_t read) { return solo->correctedTags(read); };
    return openBamSink(files, config, ref, std::move(correction));
}
} // namespace nustar::star2711b
