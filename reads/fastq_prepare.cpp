#include "fastq_prepare.hpp"
#include <algorithm>
#include <stdexcept>
#include <string_view>
#include <utility>

namespace nustar::star2711b {

namespace {
BaseCode encode(char base) {
    switch (base) {
    case 'A': case 'a': return BaseCode::A;
    case 'C': case 'c': return BaseCode::C;
    case 'G': case 'g': return BaseCode::G;
    case 'T': case 't': return BaseCode::T;
    default: return BaseCode::N;
    }
}

} // namespace

PreparedSingleRead prepareSingleFastq(std::string header, std::string sequence,
                           std::string quality, std::uint64_t ordinal,
    std::array<std::uint64_t,2> clipping,std::uint32_t file_group) {
    const auto name_end = header.find_first_of(" \t");
    auto name = header.substr(1, name_end == std::string::npos ?
                                   std::string::npos : name_end-1);
    MappingRead read{};
    read.encoded[0].reserve(sequence.size());
    read.encoded[1].reserve(sequence.size());
    read.fragment = FragmentId{file_group, ordinal};
    const auto left=std::min<std::uint64_t>(sequence.size(),clipping[0]);
    const auto right=std::min<std::uint64_t>(sequence.size()-left,clipping[1]);
    const auto length=sequence.size()-left-right;
    read.mates.push_back(PreparedMate{sequence.size(), length,
        ClipState{left,right,0,0,0,0,left,right},
        ReadRange{ReadOffset{0}, ReadOffset{length}}, false});
    for (const auto base : std::string_view(sequence).substr(left,length)) read.encoded[0].push_back(encode(base));
    for (const auto base : read.encoded[0]) {
        const auto code = static_cast<std::uint8_t>(base);
        read.encoded[1].push_back(code < 4 ? static_cast<BaseCode>(code^3) : base);
    }
    read.encoded[2].assign(read.encoded[1].rbegin(), read.encoded[1].rend());
    read.budget = ReadBudget{std::min<std::uint64_t>(10,
        std::uint64_t(length*0.3)), length, sequence.size(),
        length};
    RawMate raw{std::move(header), std::move(sequence), quality, {}};
    RawSingleRead original{read.fragment, FastxFormat::Fastq,
                           std::move(raw), 'N'};
    return PreparedSingleRead{std::move(original), std::move(name), std::move(quality),
                              SingleMappingRead::fromPrepared(std::move(read))};
}


PreparedReadPair preparePairedFastq(std::string header1, std::string sequence1,
                         std::string quality1, std::string header2,
                         std::string sequence2,
                         std::string quality2, std::uint64_t ordinal,
    std::array<std::array<std::uint64_t,2>,2> clipping,std::uint32_t file_group) {
    const auto name_end = header1.find_first_of(" \t");
    auto name = header1.substr(1, name_end == std::string::npos ?
                                   std::string::npos : name_end-1);
    MappingRead read{};
    std::array<std::string_view,2> sequences{sequence1,sequence2};
    std::array<ClipState,2> clips{};
    for(unsigned m=0;m<2;++m) {
        const auto left=std::min<std::uint64_t>(sequences[m].size(),clipping[m][0]);
        const auto right=std::min<std::uint64_t>(sequences[m].size()-left,clipping[m][1]);
        clips[m]=ClipState{left,right,0,0,0,0,left,right};
        sequences[m]=sequences[m].substr(left,sequences[m].size()-left-right);
    }
    const auto combined_length = sequences[0].size()+sequences[1].size()+1;
    read.mates.reserve(2);
    read.encoded[0].reserve(combined_length);
    read.encoded[1].reserve(combined_length);
    read.fragment = FragmentId{file_group, ordinal};
    read.mates.push_back(PreparedMate{sequence1.size(), sequences[0].size(), clips[0],
        ReadRange{ReadOffset{0}, ReadOffset{sequences[0].size()}}, false});
    for (const auto base : sequences[0]) read.encoded[0].push_back(encode(base));
    read.encoded[0].push_back(BaseCode::MateSpacer);
    for (auto it = sequences[1].rbegin(); it != sequences[1].rend(); ++it) {
        const auto code = static_cast<std::uint8_t>(encode(*it));
        read.encoded[0].push_back(code < 4 ? static_cast<BaseCode>(code^3) :
                                  BaseCode::N);
    }
    read.mates.push_back(PreparedMate{sequence2.size(), sequences[1].size(), clips[1],
        ReadRange{ReadOffset{sequences[0].size()+1},
                  ReadOffset{read.encoded[0].size()}}, true});
    for (const auto base : read.encoded[0]) {
        const auto code = static_cast<std::uint8_t>(base);
        read.encoded[1].push_back(code < 4 ? static_cast<BaseCode>(code^3) : base);
    }
    read.encoded[2].assign(read.encoded[1].rbegin(), read.encoded[1].rend());
    const auto biological = sequences[0].size()+sequences[1].size();
    read.budget = ReadBudget{std::min<std::uint64_t>(10,
        std::uint64_t(biological*0.3)), read.encoded[0].size(), sequence1.size()+sequence2.size(),
        biological};
    RawReadPair original{read.fragment, FastxFormat::Fastq,
        {RawMate{std::move(header1), std::move(sequence1), quality1, {}},
         RawMate{std::move(header2), std::move(sequence2), quality2, {}}}, 'N'};
    return PreparedReadPair{std::move(original), std::move(name), {std::move(quality1), std::move(quality2)},
                            PairedMappingRead::fromPrepared(std::move(read))};
}

} // namespace nustar::star2711b
