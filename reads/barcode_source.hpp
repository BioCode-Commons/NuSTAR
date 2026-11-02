#pragma once
#include "raw.hpp"

namespace nustar::star2711b {
// Keeps barcode FASTQ framing in the input stage; workers receive owned bytes.
class BarcodeReadSource final : public IReadSource {
public:
    BarcodeReadSource(std::unique_ptr<IReadSource> reads,std::unique_ptr<IReadSource> barcodes,
                      bool check_barcode_end)
        : reads_(std::move(reads)),barcodes_(std::move(barcodes)),check_barcode_end_(check_barcode_end) {}
    Result<std::optional<ReadBatch>> next(std::uint64_t bytes) override {
        auto result=reads_->next(bytes);
        if(auto failure=std::get_if<Failure>(&result)) return *failure;
        auto& batch=std::get<std::optional<ReadBatch>>(result);
        if(!batch && check_barcode_end_) {
            auto more=nextBarcode();
            if(auto failure=std::get_if<Failure>(&more)) return *failure;
            if(std::get<std::optional<RawMate>>(more)) return mismatch();
        } else if(batch) for(auto& read:batch->fragments) {
            auto next=nextBarcode();
            if(auto failure=std::get_if<Failure>(&next)) return *failure;
            auto& barcode=std::get<std::optional<RawMate>>(next);
            if(!barcode) return mismatch();
            read.barcode=std::move(barcode);
        }
        return result;
    }
    void cancel() noexcept override {reads_->cancel();barcodes_->cancel();}
private:
    static Failure mismatch() {return {ErrorKind::Input,Stage::Input,{},1,"RNA and barcode FASTQ record counts differ"};}
    Result<std::optional<RawMate>> nextBarcode() {
        if(!batch_ || next_==batch_->fragments.size()) {
            auto result=barcodes_->next(1ULL<<20);
            if(auto failure=std::get_if<Failure>(&result)) return *failure;
            batch_=std::get<std::optional<ReadBatch>>(std::move(result));next_=0;
            if(!batch_) return std::optional<RawMate>{};
        }
        return std::move(batch_->fragments[next_++].mates[0]);
    }
    std::unique_ptr<IReadSource> reads_,barcodes_;
    std::optional<ReadBatch> batch_;
    std::size_t next_=0;
    bool check_barcode_end_;
};
inline Result<std::unique_ptr<IReadSource>> withBarcode(std::unique_ptr<IReadSource> reads,
    const IReadSourceBackend& backend,const std::optional<std::filesystem::path>& path,
    bool check_barcode_end=true) {
    if(!path) return reads;
    InputSpec spec{Layout::SE,{},{},{}};
    spec.groups.push_back(InputGroup{{*path,{}}});
    auto barcode=backend.open(spec);
    if(auto failure=std::get_if<Failure>(&barcode)) return *failure;
    return std::unique_ptr<IReadSource>(std::make_unique<BarcodeReadSource>(std::move(reads),
        std::get<std::unique_ptr<IReadSource>>(std::move(barcode)),check_barcode_end));
}
}
