#pragma once
#include "../reads/raw.hpp"
#include "../sam/model.hpp"
#include <algorithm>
#include <string_view>
namespace nustar::star2711b {
inline void appendInputTags(ProjectedOutput& output, const RawFragment& input,
                            const std::vector<std::string>& keep, bool bam_output) {
    if (input.mates[0]->header_extra.empty() &&
        (!input.mates[1] || input.mates[1]->header_extra.empty()))
        return;
    for (auto& record : output.records) {
        const auto& text = input.mates[static_cast<unsigned>(record.mate)]->header_extra;
        std::size_t begin = 0;
        while (begin < text.size()) {
            const auto end = text.find('\t', begin);
            const auto tag = std::string_view(text).substr(
                begin, end == text.npos ? text.size() - begin : end - begin);
            if (tag.size() >= 5) {
                // STAR keeps all input tags in text SAM; its BAM converter
                // supports scalar i/A/Z/f and applies readFilesSAMattrKeep.
                const bool selected =
                    !bam_output ||
                    (keep[0] != "None" &&
                     (keep[0] == "All" ||
                      std::find(keep.begin(), keep.end(), tag.substr(0, 2)) != keep.end()));
                if (selected && (!bam_output ||
                                 std::string_view("iAZf").find(tag[3]) != std::string_view::npos))
                    record.tags.push_back({{tag[0], tag[1]}, tag[3], std::string(tag.substr(5))});
            }
            if (end == text.npos)
                break;
            begin = end + 1;
        }
    }
}
} // namespace nustar::star2711b
