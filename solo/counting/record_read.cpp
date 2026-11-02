#include "read_recorder.hpp"

namespace nustar::star2711b::solo::detail {

void ReadRecorder::record(UnsignedWord nTr, Alignment** alignOut, UnsignedWord iRead,
                          ReadAnnotations& readAnnot) {
    if (options.type == options.SoloTypes::None)
        return;
    if (options.type == options.SoloTypes::CB_samTagOut)
        return;

    if (options.readStats.yes)
        readFlagReset();

    for (std::uint32_t ii = 0; ii < options.nFeatures; ii++)
        readFeat[ii]->record(*readBar, nTr, alignOut, iRead, readAnnot);
};

} // namespace nustar::star2711b::solo::detail
