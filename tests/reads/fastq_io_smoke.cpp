#include "../../reads/raw.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>

using namespace nustar::star2711b;

namespace {

void expect(bool good, const char* message) {
    if (!good) throw std::runtime_error(message);
}

void write(const std::filesystem::path& path, const std::string& contents) {
    std::ofstream output(path, std::ios::binary);
    output << contents;
    expect(static_cast<bool>(output), "fixture write failed");
}

std::unique_ptr<IReadSource> open(const IReadSourceBackend& backend,
                                   const InputSpec& spec) {
    auto result = backend.open(spec);
    if (const auto* error = std::get_if<Failure>(&result))
        throw std::runtime_error(error->message);
    return std::get<std::unique_ptr<IReadSource>>(std::move(result));
}

std::optional<ReadBatch> next(IReadSource& source) {
    auto result = source.next(1);
    if (const auto* error = std::get_if<Failure>(&result))
        throw std::runtime_error(error->message);
    return std::get<std::optional<ReadBatch>>(std::move(result));
}

void check(const IReadSourceBackend& backend, const std::filesystem::path& dir) {
    InputSpec single{Layout::SE,
        {InputGroup{{dir/"one.fq", {}}}, InputGroup{{dir/"two.fq", {}}}}, {}, {}};
    auto source = open(backend, single);
    for (std::uint64_t ordinal = 0; ordinal < 3; ++ordinal) {
        auto batch = next(*source);
        expect(batch && batch->id.value == ordinal && batch->fragments.size() == 1,
               "single-end batch order failed");
        expect(batch->fragments.front().id.ordinal == ordinal,
               "single-end ordinal failed");
        expect(batch->fragments.front().mates[0]->quality->size() == 2,
               "single-end quality failed");
    }
    expect(!next(*source), "single-end EOF failed");
    expect(!next(*source), "repeated EOF failed");

    InputSpec pair{Layout::PE, {InputGroup{{dir/"one.fq", dir/"mates.fq"}}}, {}, {}};
    auto paired = open(backend, pair);
    auto first = next(*paired);
    expect(first && first->fragments.front().mates[1].has_value(),
           "paired fragment missing mate");
    auto incomplete = paired->next(1);
    expect(std::holds_alternative<Failure>(incomplete),
           "unmatched mate EOF was accepted");

    InputSpec empty{Layout::SE, {InputGroup{{dir/"empty.fq", {}}}}, {}, {}};
    auto empty_source = open(backend, empty);
    expect(!next(*empty_source), "empty FASTQ did not reach EOF");

    InputSpec separators{Layout::SE, {InputGroup{{dir/"separators.fq", {}}}}, {}, {}};
    auto permissive = open(backend, separators);
    for (unsigned i = 0; i < 2; ++i) {
        auto batch = next(*permissive);
        expect(batch && batch->fragments.size() == 1 &&
               batch->fragments[0].mates[0]->sequence == "AC" &&
               *batch->fragments[0].mates[0]->quality == "!!",
               "STAR-compatible third-line consumption failed");
    }
    expect(!next(*permissive), "third-line fixture has extra records");
    for (const auto* name : {"bad-header.fq", "bad-quality.fq", "truncated.fq"}) {
        InputSpec invalid{Layout::SE, {InputGroup{{dir/name, {}}}}, {}, {}};
        auto rejected = open(backend, invalid);
        expect(std::holds_alternative<Failure>(rejected->next(1)),
               "required FASTQ framing or quality validation was removed");
    }
}

} // namespace

int main(int argc, char** argv) {
    try {
        expect(argc == 2, "expected fixture directory");
        const std::filesystem::path dir(argv[1]);
        std::filesystem::create_directories(dir);
        write(dir/"one.fq", "@r1\r\nAC\r\n+\r\n!!\r\n@r2\nGT\n+\n##");
        write(dir/"two.fq", "@r3\nTA\n+\n$$\n");
        write(dir/"mates.fq", "@r1\nGT\n+\n!!\n");
        write(dir/"empty.fq", "");
        write(dir/"separators.fq", "@r1\nAC\nignored\n!!\n@r2\nAC\n\n!!\n");
        write(dir/"bad-header.fq", "r1\nAC\n+\n!!\n");
        write(dir/"bad-quality.fq", "@r1\nAC\n+\n!\n");
        write(dir/"truncated.fq", "@r1\nAC\n+\n");
        check(memoryMappedFastqBackend(), dir);
        check(bufferedFastqBackend(), dir);
        // Cross several mmap release boundaries. Keep an earlier batch alive
        // while reading later data to verify that returned strings own bytes.
        const std::string sequence(151, 'A'), quality(151, 'I');
        {
            std::ofstream output(dir/"large.fq", std::ios::binary);
            for (unsigned n = 0; n < 120000; ++n)
                output << "@large" << n << '\n' << sequence << "\n+\n" << quality << '\n';
            expect(static_cast<bool>(output), "large fixture write failed");
        }
        InputSpec spec{Layout::SE, {InputGroup{{dir/"large.fq", {}}}}, {}, {}};
        auto mapped = open(memoryMappedFastqBackend(), spec);
        auto buffered = open(bufferedFastqBackend(), spec);
        auto retained = next(*mapped);
        auto reference = next(*buffered);
        expect(retained && reference, "large fixture first read missing");
        for (unsigned n = 1; n < 120000; ++n) {
            auto a = next(*mapped), b = next(*buffered);
            expect(a && b, "large fixture truncated");
            const auto& left = *a->fragments[0].mates[0];
            const auto& right = *b->fragments[0].mates[0];
            expect(left.header == right.header && left.sequence == right.sequence &&
                   left.quality == right.quality, "mmap release changed FASTQ contents");
        }
        expect(!next(*mapped) && !next(*buffered), "large fixture EOF failed");
        expect(retained->fragments[0].mates[0]->sequence == sequence &&
               retained->fragments[0].mates[0]->quality == quality,
               "retained batch lost ownership after mmap release");
        std::filesystem::remove(dir/"large.fq");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
