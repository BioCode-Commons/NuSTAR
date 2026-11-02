#include "encode.hpp"

#include <algorithm>
#include <charconv>
#include <sstream>

namespace nustar::star2711b {
namespace {

std::vector<std::byte> bytes(std::string_view text) {
    if (text.empty()) return {};
    const auto* first = reinterpret_cast<const std::byte*>(text.data());
    return {first, first+text.size()};
}

template<class Number> void appendNumber(std::string& output, Number value) {
    char buffer[32];
    const auto result = std::to_chars(buffer, buffer+sizeof(buffer), value);
    output.append(buffer, result.ptr);
}

void appendRecord(std::string& output, const SamRecord& record) {
    output += record.qname;
    output += '\t'; appendNumber(output, record.flag);
    output += '\t'; output += record.rname;
    output += '\t'; appendNumber(output, record.position_1based);
    output += '\t'; appendNumber(output, unsigned(record.mapq));
    output += '\t';
    if (record.cigar.empty()) output += '*';
    else for (const auto& operation : record.cigar) {
        appendNumber(output, operation.length);
        output += operation.operation;
    }
    output += '\t'; output += record.rnext;
    output += '\t'; appendNumber(output, record.mate_position_1based);
    output += '\t'; appendNumber(output, record.template_length);
    output += '\t'; output += record.sequence;
    output += '\t'; output += record.quality;
    for (const auto& tag : record.tags) {
        output += '\t'; output += tag.tag[0]; output += tag.tag[1];
        output += ':'; output += tag.type; output += ':';
        output += tag.encoded_value;
    }
    output += '\n';
}

} // namespace

Result<EncodedFragment> encodeSam(const FragmentId& fragment,
                                  const ProjectedOutput& projected) {
    std::string output;
    std::size_t estimate = 0;
    for (const auto& record : projected.records)
        estimate += record.sequence.size()+record.quality.size()+record.qname.size()+128;
    output.reserve(estimate);
    for (const auto& record : projected.records) {
        appendRecord(output, record);
    }
    EncodedFragment encoded{};
    encoded.fragment = fragment;
    encoded.sam = bytes(output);
    for (const auto& record : projected.unmapped)
        encoded.unmapped_fastx[static_cast<unsigned>(record.mate)] =
            bytes(record.encoded_record);
    return encoded;
}

Result<std::vector<std::byte>> encodeSamHeader(const ReferenceMetadata& reference,
    const SamConfig& config, const BaselineIdentity& baseline) {
    if (baseline.release != "2.7.11b")
        return Failure{ErrorKind::Config, Stage::Project, {}, 1,
                       "SAM header requires a STAR 2.7.11b baseline"};
    std::ostringstream output;
    output << (config.header_hd.empty() ? "@HD\tVN:1.4" : config.header_hd)
           << '\n';
    for (const auto& chromosome : reference.chromosomes)
        output << "@SQ\tSN:" << chromosome.name << "\tLN:"
               << chromosome.length << '\n';
    for (const auto& line : config.extra_reference_lines)
        output << line << '\n';
    if (!config.header_pg.empty() && config.header_pg != "-")
        output << config.header_pg << '\n';
    output << "@PG\tID:STAR\tPN:STAR\tVN:" << baseline.release
           << "\tCL:" << config.command_line_full << '\n';
    for (const auto& comment : config.header_comments)
        output << comment << '\n';
    for (const auto& read_group : config.read_group_lines)
        output << "@RG\t" << read_group << '\n';
    output << "@CO\tuser command line: " << config.command_line_user << '\n';
    output << config.header_extra;
    return bytes(output.str());
}


} // namespace nustar::star2711b
