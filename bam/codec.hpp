#pragma once
#include "../reference/view.hpp"
#include "../sam/model.hpp"
#include <cstddef>
#include <span>
#include <unordered_map>
namespace nustar::star2711b::bam {
using Bytes = std::vector<std::byte>;
void append32(Bytes&, std::uint32_t);
std::uint32_t read32(const std::byte*);
void appendStringTag(Bytes&, const char*, const std::string&);
class Encoder {
    std::unordered_map<std::string, std::int32_t> references_;

  public:
    explicit Encoder(const ReferenceMetadata&);
    void append(Bytes&, const SamRecord&) const;
};
} // namespace nustar::star2711b::bam
