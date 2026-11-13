#pragma once

#include "raw.hpp"
#include <algorithm>

namespace nustar::star2711b {

class LimitedReadSource final : public IReadSource {
public:
    LimitedReadSource(std::unique_ptr<IReadSource> source, std::uint64_t limit)
        : source_(std::move(source)), remaining_(limit), limit_(limit) {}

    Status prepareParallel(std::uint64_t bytes, std::uint32_t threads,
                           std::uint64_t limit = UINT64_MAX) override {
        auto result = source_->prepareParallel(bytes, threads, std::min(limit, limit_));
        if (std::holds_alternative<std::monostate>(result))
            parallel_ = source_->concurrentReads();
        return result;
    }
    bool concurrentReads() const noexcept override { return parallel_; }

    Result<std::optional<ReadBatch>> next(std::uint64_t bytes) override {
        if (parallel_) {
            // The parallel source applies the limit when claiming record batches.
            return source_->next(bytes);
        }
        if (!remaining_) return std::optional<ReadBatch>{};
        auto result = source_->next(bytes);
        if (auto* batch = std::get_if<std::optional<ReadBatch>>(&result); batch && *batch) {
            auto& reads = (**batch).fragments;
            if (reads.size() > remaining_) reads.resize(remaining_);
            remaining_ -= reads.size();
        }
        return result;
    }
    void cancel() noexcept override { source_->cancel(); }
private:
    std::unique_ptr<IReadSource> source_;
    std::uint64_t remaining_;
    const std::uint64_t limit_;
    bool parallel_ = false;
};

inline std::unique_ptr<IReadSource> limitReads(std::unique_ptr<IReadSource> source,
                                             std::uint64_t limit) {
    if (limit == UINT64_MAX) return source;
    return std::make_unique<LimitedReadSource>(std::move(source), limit);
}

} // namespace nustar::star2711b
