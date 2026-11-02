#pragma once

#include "../annotation/read_annotations.hpp"
#include "../barcode/read_barcode.hpp"
#include "feature_recorder.hpp"

namespace nustar::star2711b::solo::detail {

class ReadRecorder {
public:
    std::unique_ptr<ReadBarcode> readBar;
    std::vector<std::unique_ptr<FeatureRecorder>> readFeat;

    ReadRecorder(Context& context_in, std::int32_t worker_id_in);
    ~ReadRecorder() = default;
    ReadRecorder(const ReadRecorder&) = delete;
    ReadRecorder& operator=(const ReadRecorder&) = delete;
    void readFlagReset();
    void record(UnsignedWord nTr, Alignment** alignOut, UnsignedWord iRead,
                ReadAnnotations& readAnnot);

private:
    const std::int32_t worker_id;
    Context& context_;
    Options& options;
};

} // namespace nustar::star2711b::solo::detail
