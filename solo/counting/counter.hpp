#pragma once

#include "../annotation/index.hpp"
#include "feature_counter.hpp"
#include "worker_state.hpp"
#include "../support/types.hpp"
#include <fstream>

namespace nustar::star2711b::solo::detail {

class Counter {
private:
    WorkerState** workers;
    Context& context_;
    AnnotationIndex& annotation_;

public:
    Options& options;
    FeatureCounter** features_ = nullptr;

    ReadBarcode* barcode_sum_ = nullptr;

    Counter(WorkerState** workers, Context& context_in, AnnotationIndex& annotation_in);

    ~Counter();
    Counter(const Counter&) = delete;
    Counter& operator=(const Counter&) = delete;

    void processAndOutput();
};

} // namespace nustar::star2711b::solo::detail
