#include "read_recorder.hpp"

namespace nustar::star2711b::solo::detail {

ReadRecorder::ReadRecorder(Context& context_in, std::int32_t worker_id_in)
    : worker_id(worker_id_in), context_(context_in), options(context_.options) {
    readBar = std::make_unique<ReadBarcode>(context_);

    if (options.type == 0)
        return;
    if (options.type == options.SoloTypes::CB_samTagOut)
        return;

    readFeat.resize(options.nFeatures);

    for (std::uint32_t ii = 0; ii < options.nFeatures; ii++)
        readFeat[ii] = std::make_unique<FeatureRecorder>(options.features[ii], context_, worker_id);
};

void ReadRecorder::readFlagReset() {
    for (std::uint32_t ii = 0; ii < options.nFeatures; ii++)
        readFeat[ii]->readFlag.flag = 0;
};

} // namespace nustar::star2711b::solo::detail
