#pragma once

#include "raw.hpp"

namespace nustar::star2711b {

class LimitedReadSource final : public IReadSource {
public:
    LimitedReadSource(std::unique_ptr<IReadSource> source, std::uint64_t limit)
        : source_(std::move(source)), remaining_(limit) {}

    Result<std::optional<ReadBatch>> next(std::uint64_t bytes) override {
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
};

inline std::unique_ptr<IReadSource> limitReads(std::unique_ptr<IReadSource> source,
                                             std::uint64_t limit) {
    if (limit == UINT64_MAX) return source;
    return std::make_unique<LimitedReadSource>(std::move(source), limit);
}

} // namespace nustar::star2711b
