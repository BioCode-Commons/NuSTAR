#include "counter.hpp"
#include "../support/streams.hpp"
#include "../support/time.hpp"

namespace nustar::star2711b::solo::detail {

Counter::Counter(WorkerState** workers_in, Context& context_in, AnnotationIndex& annotation_in)
    : workers(workers_in), context_(context_in), annotation_(annotation_in),
      options(context_.options) {
    if (options.type == 0)
        return;

    barcode_sum_ = new ReadBarcode(context_);

    if (options.type == options.SoloTypes::CB_samTagOut)
        return;

    features_ = new FeatureCounter*[options.nFeatures];
    for (std::uint32_t ii = 0; ii < options.nFeatures; ii++)
        features_[ii] = new FeatureCounter(context_, workers, annotation_, options.features[ii],
                                           barcode_sum_, features_);
};

void Counter::processAndOutput() {
    if (options.type == 0)
        return;

    {                          // collect barcode statistics
        if (options.cbWLyes) { // now we can define WL and counts
            for (int ii = 0; ii < context_.runThreadN; ii++) {
                barcode_sum_->addCounts(*workers[ii]->read->readBar);
                barcode_sum_->addStats(*workers[ii]->read->readBar);
                workers[ii]->read->readBar.reset(); // no longer needed
            };
        };

        std::ofstream* statsStream =
            &ofstrOpen(context_.outFileNamePrefix + options.outFileNames[0] + "Barcodes.stats",
                       sourceLocation(), context_);
        barcode_sum_->statsOut(*statsStream);
        statsStream->close();

        // pseudocounts
        if (options.CBmatchWL.mm1_multi_pc) {
            for (std::uint32_t ii = 0; ii < options.cbWLsize; ii++) {
                barcode_sum_->cbReadCountExact[ii]++; // add one to exact counts
            };
        };
    };

    if (options.type == options.SoloTypes::CB_samTagOut)
        return;

    { // process all features
        *context_.inOut->logStdOut << timeMonthDayTime() << " ..... started Solo counting\n"
                                   << std::flush;
        context_.inOut->logMain << timeMonthDayTime() << " ..... started Solo counting\n"
                                << std::flush;

        for (std::uint32_t ii = 0; ii < options.nFeatures; ii++) {
            features_[ii]->processRecords();
            //};
        };

        *context_.inOut->logStdOut << timeMonthDayTime() << " ..... finished Solo counting\n"
                                   << std::flush;
        context_.inOut->logMain << timeMonthDayTime() << " ..... finished Solo counting\n"
                                << std::flush;
    };
};

Counter::~Counter() {
    if (features_) {
        for (std::uint32_t i = 0; i < options.nFeatures; ++i)
            delete features_[i];
        delete[] features_;
    }
    delete barcode_sum_;
}

} // namespace nustar::star2711b::solo::detail
