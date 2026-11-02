#pragma once

#include "../annotation/index.hpp"
#include "read_recorder.hpp"
#include "record_types.hpp"
#include "worker_state.hpp"
#include "../filtering/cell_statistics.hpp"
#include "../support/types.hpp"
#include <fstream>
#include <unordered_map>
#include <unordered_set>

namespace nustar::star2711b::solo::detail {

class FeatureCounter {
private:
    Context& context_;
    WorkerState** workers;
    AnnotationIndex& annotation_;

    FeatureCounter** features_all_;

    static const std::uint32_t umiArrayStride = 3;
    enum { rguG, rguU, rguR };
    std::uint32_t rguStride;

public:
    Options& options;

    FeatureRecorder *readFeatSum = nullptr, **readFeatAll = nullptr;
    ReadBarcode* barcode_sum_;

    const std::int32_t featureType;

    UnsignedWord nReadsMapped, nReadsInput; // total number of mapped reads
    std::uint32_t nCB;
    std::uint32_t featuresNumber; // number of features (i.e. genes, SJs, etc)

    std::uint32_t* rGeneUMI; // mapped reads sorted by CB
    std::uint32_t* rCBn;     // number of reads for detected CBs in the whitelist
    std::uint32_t** rCBp;    // array of pointers to each CB sub-array

    std::vector<std::uint32_t> indCB;   // index of detected CBs in the whitelist
    std::vector<std::uint32_t> indCBwl; // reverse of indCB: index of WL CBs in detected CB list
    std::vector<std::uint32_t> nUMIperCB,
        nUMIperCBsorted; // number of UMIs per CB, and the same sorted (descendant)
    std::vector<std::uint32_t> nGenePerCB, nGenePerCBmulti; // number of genes (with >0 UMIs) per CB
    std::vector<std::uint32_t>
        nReadPerCB; // number of reads per CB. With multimappers: all aligns per CB
    std::vector<std::uint32_t> nReadPerCBunique,
        nReadPerCBtotal; // number of unique and multiple reads per CB
    std::uint32_t nReadPerCBmax;

    std::vector<double> nUMIperCBmulti;

    std::vector<std::uint32_t> countCellGeneUMI; // sparsified matrix for the counts, each entry is:
                                                 // geneID count1 count2 ... countNcounts
    std::vector<std::uint32_t> countCellGeneUMIindex; // index of CBs in the count matrix
    std::uint32_t countMatStride; // number of counts per entry in the count matrix

    struct {
        std::vector<double> m;
        std::vector<std::uint32_t> i;
        std::uint32_t s;
    } countMatMult;

    std::vector<std::unordered_map<std::uint32_t, std::unordered_set<UnsignedWord>>>
        cbFeatureUMImap; // for SmartSeq counting

    std::string outputPrefix, outputPrefixFiltered;

    CellStatistics filteredCells;

    std::array<std::vector<UnsignedWord>, 2> sjAll;

    std::vector<ReadCorrection> readInfo; // corrected CB/UMI information for each read
    ReadFlags readFlagCounts;

    std::vector<std::uint32_t> redistrFilesCBindex,
        redistrFilesCBfirst; // redistr file for each CB, CB boundaries in redistributed files
    std::vector<UnsignedWord> redistrFilesNreads; // number of reads in each file
    std::vector<std::fstream*> redistrFilesStreams;

    FeatureCounter(Context& context_in, WorkerState** workers, AnnotationIndex& annotation_in,
                   std::int32_t feTy, ReadBarcode* readBarSumIn, FeatureCounter** features_all_);
    ~FeatureCounter();
    FeatureCounter(const FeatureCounter&) = delete;
    FeatureCounter& operator=(const FeatureCounter&) = delete;
    void clearLarge(); // clear large vectors
    void processRecords();
    void sumThreads();
    void countSmartSeq();
    void countCBgeneUMI();
    void countVelocyto();
    void quantTranscript();

    void collapseUMI(std::uint32_t iCB, std::uint32_t* umiArray);
    void collapseUMI_CR(std::uint32_t iCB, std::uint32_t* umiArray);
    void deduplicateAllCells();
    void deduplicateCell(std::uint32_t iCB, std::vector<std::uint32_t>& umiArray,
                         std::vector<std::uint32_t>& gID, std::vector<std::uint32_t>& gReadS);

    std::uint32_t correctUmisCellRanger(const std::uint32_t nU0, UmiCode* umiArr,
                                        const bool readInfoRec, const bool nUMIyes,
                                        std::unordered_map<UmiCode, UmiCode>& umiCorr);
    std::uint32_t correctUmisDirectional(const std::uint32_t nU0, UmiCode* umiArr,
                                         const bool readInfoRec, const bool nUMIyes,
                                         std::unordered_map<UmiCode, UmiCode>& umiCorr,
                                         const std::int32_t dirCountAdd);
    std::uint32_t correctUmisGraph(const std::uint32_t nU0, UmiCode* umiArr, const bool readInfoRec,
                                   const bool nUMIyes,
                                   std::unordered_map<UmiCode, UmiCode>& umiCorr);

    void outputResults(bool cellFilterYes, std::string outputPrefixMat);

    void statsOutput();
    void redistributeReadsByCB();

    void cellFiltering();
    void emptyDrops_CR();
    void loadRawMatrix();
};

} // namespace nustar::star2711b::solo::detail
