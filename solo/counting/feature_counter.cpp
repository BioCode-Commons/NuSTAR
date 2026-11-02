#include "feature_counter.hpp"
#include "../support/streams.hpp"

namespace nustar::star2711b::solo::detail {

FeatureCounter::FeatureCounter(Context& context_in, WorkerState** workers,
                               AnnotationIndex& annotation_in, std::int32_t feTy,
                               ReadBarcode* readBarSumIn, FeatureCounter** features_all_)
    : context_(context_in), workers(workers), annotation_(annotation_in), featureType(feTy),
      features_all_(features_all_), options(context_.options), barcode_sum_(readBarSumIn) {
    if (featureType >=
        0) { // otherwise we do not need these arrays - e.g. with --runMode soloCellFiltering
        readFeatSum = new FeatureRecorder(featureType, context_, -1);
        readFeatAll = new FeatureRecorder*[context_.runThreadN];
    };

    // number of features
    switch (featureType) {
    case FeatureTypes::Gene:
    case FeatureTypes::GeneFull:
    case FeatureTypes::GeneFull_Ex50pAS:
    case FeatureTypes::GeneFull_ExonOverIntron:
    case FeatureTypes::Velocyto:
        featuresNumber = annotation_.nGe;
        break;
    case FeatureTypes::SJ:
        featuresNumber = context_.sjAll[0].size();
        break;
    default:
        featuresNumber = -1; // undefined
    };
};

void FeatureCounter::clearLarge() {
    cbFeatureUMImap.clear();
    cbFeatureUMImap.shrink_to_fit();
    countCellGeneUMI.clear();
    countCellGeneUMI.shrink_to_fit();
    countCellGeneUMIindex.clear();
    countCellGeneUMIindex.shrink_to_fit();
    countMatMult.i.clear();
    countMatMult.i.shrink_to_fit();
    countMatMult.m.clear();
    countMatMult.m.shrink_to_fit();
    // indCB.clear(); //needed for Velocyto
    // indCB.shrink_to_fit();
    indCBwl.clear();
    indCBwl.shrink_to_fit();
    nGenePerCB.clear();
    nGenePerCB.shrink_to_fit();
    nGenePerCBmulti.clear();
    nGenePerCBmulti.shrink_to_fit();
    nReadPerCB.clear();
    nReadPerCB.shrink_to_fit();
    nReadPerCBtotal.clear();
    nReadPerCBtotal.shrink_to_fit();
    nReadPerCBunique.clear();
    nReadPerCBunique.shrink_to_fit();
    nUMIperCB.clear();
    nUMIperCB.shrink_to_fit();
    nUMIperCBmulti.clear();
    nUMIperCBmulti.shrink_to_fit();
    nUMIperCBsorted.clear();
    nUMIperCBsorted.shrink_to_fit();
    sjAll[0].clear();
    sjAll[0].shrink_to_fit();
    sjAll[1].clear();
    sjAll[1].shrink_to_fit();
};

FeatureCounter::~FeatureCounter() {
    delete readFeatSum;
    delete[] readFeatAll;
}

} // namespace nustar::star2711b::solo::detail
