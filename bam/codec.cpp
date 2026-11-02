#include "codec.hpp"
#include <algorithm>
#include <bit>
#include <charconv>
#include <stdexcept>
#include <string_view>
namespace nustar::star2711b::bam {
void append32(Bytes& out, std::uint32_t n) {
    for (unsigned i = 0; i < 4; ++i)
        out.push_back(std::byte(n >> (i * 8)));
}
std::uint32_t read32(const std::byte* p) {
    return std::to_integer<std::uint32_t>(p[0]) | std::to_integer<std::uint32_t>(p[1]) << 8 |
           std::to_integer<std::uint32_t>(p[2]) << 16 | std::to_integer<std::uint32_t>(p[3]) << 24;
}
void appendStringTag(Bytes& out, const char* tag, const std::string& value) {
    out.push_back(std::byte(tag[0]));
    out.push_back(std::byte(tag[1]));
    out.push_back(std::byte('Z'));
    const auto* p = reinterpret_cast<const std::byte*>(value.c_str());
    out.insert(out.end(), p, p + value.size() + 1);
}
Encoder::Encoder(const ReferenceMetadata& ref) {
    for (std::size_t i = 0; i < ref.chromosomes.size(); ++i)
        references_.emplace(ref.chromosomes[i].name, i);
}
namespace {
unsigned base(char c) {
    switch (c) {
    case 'A':
    case 'a':
        return 1;
    case 'C':
    case 'c':
        return 2;
    case 'G':
    case 'g':
        return 4;
    case 'T':
    case 't':
        return 8;
    case '=':
        return 0;
    case 'M':
        return 3;
    case 'R':
        return 5;
    case 'S':
        return 6;
    case 'V':
        return 7;
    case 'W':
        return 9;
    case 'Y':
        return 10;
    case 'H':
        return 11;
    case 'K':
        return 12;
    case 'D':
        return 13;
    case 'B':
        return 14;
    default:
        return 15;
    }
}
unsigned bin(std::int32_t start, std::int32_t end) {
    --end;
    if (start >> 14 == end >> 14)
        return 4681 + (start >> 14);
    if (start >> 17 == end >> 17)
        return 585 + (start >> 17);
    if (start >> 20 == end >> 20)
        return 73 + (start >> 20);
    if (start >> 23 == end >> 23)
        return 9 + (start >> 23);
    if (start >> 26 == end >> 26)
        return 1 + (start >> 26);
    return 0;
}
} // namespace
void Encoder::append(Bytes& out, const SamRecord& r) const {
    const auto begin = out.size();
    const auto tid = r.rname == "*" ? -1 : references_.at(r.rname);
    const auto mtid = r.rnext == "*" ? -1 : r.rnext == "=" ? tid : references_.at(r.rnext);
    const auto pos = std::int32_t(r.position_1based) - 1;
    const auto length = r.sequence == "*" ? 0 : r.sequence.size();
    if (r.qname.size() > 254 || r.cigar.size() > 65535)
        throw std::runtime_error("BAM QNAME/CIGAR exceeds format limit");
    auto end = pos;
    for (const auto& c : r.cigar)
        if (std::string_view("MDN=X").find(c.operation) != std::string_view::npos)
            end += c.length;
    append32(out, 0);
    append32(out, tid);
    append32(out, pos);
    append32(out, (bin(pos, std::max(pos + 1, end)) << 16) | (unsigned(r.mapq) << 8) |
                      (r.qname.size() + 1));
    append32(out, (unsigned(r.flag) << 16) | r.cigar.size());
    append32(out, length);
    append32(out, mtid);
    append32(out, std::int32_t(r.mate_position_1based) - 1);
    append32(out, r.template_length);
    const auto* q = reinterpret_cast<const std::byte*>(r.qname.c_str());
    out.insert(out.end(), q, q + r.qname.size() + 1);
    for (const auto& c : r.cigar)
        append32(out, (c.length << 4) | std::string_view("MIDNSHP=XB").find(c.operation));
    for (std::size_t i = 0; i < length; i += 2)
        out.push_back(
            std::byte((base(r.sequence[i]) << 4) | (i + 1 < length ? base(r.sequence[i + 1]) : 0)));
    for (std::size_t i = 0; i < length; ++i)
        out.push_back(r.quality == "*" ? std::byte(255) : std::byte(r.quality[i] - 33));
    for (const auto& t : r.tags) {
        if (t.type == 'Z') {
            appendStringTag(out, t.tag.data(), t.encoded_value);
            continue;
        }
        out.push_back(std::byte(t.tag[0]));
        out.push_back(std::byte(t.tag[1]));
        out.push_back(std::byte(t.type));
        if (t.type == 'A')
            out.push_back(std::byte(t.encoded_value[0]));
        else if (t.type == 'i') {
            std::int64_t n = 0;
            std::from_chars(t.encoded_value.data(), t.encoded_value.data() + t.encoded_value.size(),
                            n);
            unsigned width = 4;
            char type = n < 0 ? 'i' : 'I';
            if (n >= 0 && n <= 255) {
                type = 'C';
                width = 1;
            } else if (n >= 0 && n <= 65535) {
                type = 'S';
                width = 2;
            } else if (n >= -127 && n < 0) {
                type = 'c';
                width = 1;
            } else if (n >= -32767 && n < 0) {
                type = 's';
                width = 2;
            }
            out.back() = std::byte(type);
            for (unsigned i = 0; i < width; ++i)
                out.push_back(std::byte(std::uint64_t(n) >> (8 * i)));
        } else if (t.type == 'f')
            append32(out, std::bit_cast<std::uint32_t>(std::stof(t.encoded_value)));
        else if (t.type == 'B') {
            const char subtype = t.encoded_value[0];
            out.push_back(std::byte(subtype));
            append32(out, std::count(t.encoded_value.begin(), t.encoded_value.end(), ','));
            std::size_t offset = 2;
            while (offset < t.encoded_value.size()) {
                const auto comma = t.encoded_value.find(',', offset);
                const auto n = std::stoll(t.encoded_value.substr(offset, comma - offset));
                if (subtype == 'c' || subtype == 'C')
                    out.push_back(std::byte(n));
                else
                    append32(out, n);
                if (comma == std::string::npos)
                    break;
                offset = comma + 1;
            }
        } else
            throw std::runtime_error("unsupported BAM auxiliary type");
    }
    const auto n = out.size() - begin - 4;
    for (unsigned i = 0; i < 4; ++i)
        out[begin + i] = std::byte(n >> (8 * i));
}
} // namespace nustar::star2711b::bam
