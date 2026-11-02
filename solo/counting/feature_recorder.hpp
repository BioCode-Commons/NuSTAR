#pragma once

#include "../annotation/alignment.hpp"
#include "../annotation/read_annotations.hpp"
#include "../barcode/read_barcode.hpp"
#include "../configuration/context.hpp"
#include "read_statistics.hpp"
#include "record_types.hpp"
#include "../support/types.hpp"
#include <map>
#include <set>

namespace nustar::star2711b::solo::detail {

class FeatureCounter;

class FeatureRecorder {
public:
    std::uint32_t homoPolymer[4]; // homopolymer constants

    std::vector<std::uint32_t> cbReadCount;
    std::map<BarcodeCode, std::uint32_t> cbReadCountMap;

    std::vector<std::uint32_t> transcriptDistCount;

    bool readInfoYes, readIndexYes;

    std::fstream* streamReads;

    std::string cbSeq, umiSeq, cbQual, umiQual;

    ReadFlags readFlag;

    ReadStatistics stats;

    FeatureRecorder(std::int32_t feTy, Context& context_in, int worker_id);
    void record(ReadBarcode& soloBar, UnsignedWord nTr, Alignment** alignOut, UnsignedWord iRead,
                ReadAnnotations& readAnnot);
    void addCounts(const FeatureRecorder& soloCBin);
    void addStats(const FeatureRecorder& soloCBin);
    void statsOut(std::ofstream& streamOut);
    void inputRecords(std::uint32_t** cbP, std::uint32_t cbPstride,
                      std::vector<std::uint32_t>& cbReadCountTotal,
                      std::vector<ReadCorrection>& readInfo, ReadFlags& readFlagCounts,
                      std::vector<std::uint32_t>& nReadPerCBunique1,
                      std::vector<std::uint32_t>& nReadPerCBmulti1);

private:
    const std::int32_t featureType;

    Context& context_;
    Options& options;
};

} // namespace nustar::star2711b::solo::detail
