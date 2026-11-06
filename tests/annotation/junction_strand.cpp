#include "../../annotation/junction_database.hpp"

#include <algorithm>
#include <array>
#include <iostream>
#include <stdexcept>

using namespace nustar::star2711b;

namespace {
struct Sequence final : IReferenceSequence {
    std::vector<BaseCode> bases = std::vector<BaseCode>(1024, BaseCode::N);
    Result<SequenceWindow> pin(GenomeRange range) const override {
        return SequenceWindow{range, std::span(bases).subspan(
            range.begin.value, range.end.value-range.begin.value), {}};
    }
};
struct Input { Strand strand; std::uint8_t priority = 20; };

void check(std::string_view motif, std::vector<Input> input, Strand expected,
           std::size_t winner, bool separate_groups = true) {
    Sequence sequence;
    const auto base = [](char c) {
        return c == 'A' ? BaseCode::A : c == 'C' ? BaseCode::C :
               c == 'G' ? BaseCode::G : c == 'T' ? BaseCode::T : BaseCode::N;
    };
    sequence.bases[400] = base(motif[0]); sequence.bases[401] = base(motif[1]);
    sequence.bases[498] = base(motif[2]); sequence.bases[499] = base(motif[3]);
    ReferenceMetadata metadata{{{{0}, "chr1", {0}, 1024}}, {1024}, "test"};
    ReferenceView view{{1}, "test", &metadata, &sequence, {}};
    std::vector<JunctionProposal> proposals;
    const auto add = [&](std::uint64_t start, Strand strand, std::uint8_t priority) {
        proposals.push_back({{0}, {{start}, {start+100}}, strand,
            JunctionOrigin::JunctionFile, priority, proposals.size(), {}});
    };
    for (const auto value : input) add(400, value.strand, value.priority);
    // Separate the ends of STAR's strand groups. A lone duplicate at their
    // boundary is collapsed in the first pass, before opposite-strand rules.
    if (separate_groups)
        for (const auto strand : {Strand::Plus, Strand::Minus, Strand::Unknown}) {
            add(100, strand, 20); add(700, strand, 20);
        }
    AnnotationConfig config{};
    config.sjdbOverhang = 20; config.limitSjdbInsertNsj = 100;
    auto result = GenomeJunctionPreparer{}.prepare(proposals, view, config);
    if (auto* error = std::get_if<Failure>(&result)) throw std::runtime_error(error->message);
    const auto& database = std::get<PreparedJunctionDatabase>(result);
    for (std::size_t i = 0; i < database.entries.size(); ++i) {
        const auto& entry = database.entries[i];
        if (entry.reported_intron.begin.value != 400) continue;
        if (entry.strand != expected || database.winning_proposal_indices[i] != winner)
            throw std::runtime_error("unexpected junction strand/winner for " + std::string(motif));
        return;
    }
    throw std::runtime_error("target junction missing");
}
}

int main() {
    try {
        const auto p = Strand::Plus, m = Strand::Minus, u = Strand::Unknown;
        for (const auto motif : {"GTAG", "GCAG", "ATAC", "CTAC", "CTGC", "GTAT"}) {
            const bool reverse = std::string_view(motif) == "CTAC" ||
                std::string_view(motif) == "CTGC" || std::string_view(motif) == "GTAT";
            const auto correct = reverse ? m : p;
            check(motif, {{p}, {m}}, correct, reverse ? 1 : 0);
            check(motif, {{m}, {p}}, correct, reverse ? 0 : 1);
            check(motif, {{u}, {m}, {p}}, correct, reverse ? 1 : 2);
            check(motif, {{u}}, correct, 0);
        }
        check("CTAC", {{p}}, p, 0); // An explicit lone annotation is not corrected.
        check("CTAC", {{p,30}, {m,20}}, p, 0); // Priority precedes motif strand.
        check("GTAG", {{p,10}, {m,30}}, m, 1);
        check("CTAC", {{p,30}, {m,10}, {u,20}}, m, 2); // STAR's previous sorted priority.
        check("CTAC", {{p,20}, {p,30}, {m,20}}, p, 1);
        check("NNNN", {{p}, {m}}, u, 0);
        check("NNNN", {{m}, {p}}, u, 1); // Strand order, not input order.
        check("NNNN", {{p}, {m}, {u}}, u, 0);
        check("NNNN", {{p}, {u}}, p, 0);
        check("NNNN", {{m}, {u}}, m, 0);
        check("NNNN", {{p,20}, {m,30}}, m, 1);
        check("NNNN", {{p,30}, {m,20}}, p, 0);
        check("NNNN", {{u,30}, {p,20}}, u, 0);
        check("CTAC", {{p}, {m}}, p, 0, false);
        check("CTAC", {{m}, {p}}, p, 1, false);
        std::cout << "Passed 39 junction strand/priority cases\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n'; return 1;
    }
}
