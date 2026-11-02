#pragma once
#include "annotation/index.hpp"
#include "counting/worker_state.hpp"
#include "session.hpp"

namespace nustar::star2711b::solo {
struct Worker::Impl {
    detail::Context& parameters;
    detail::AnnotationIndex& transcriptome;
    detail::WorkerState native{};
    StatsDelta statistics{};
    std::vector<detail::Alignment> transcripts;
    std::vector<detail::Alignment*> alignments;
    detail::ReadAnnotations annotations;
    Impl(detail::Context& p, detail::AnnotationIndex& t, unsigned id)
        : parameters(p), transcriptome(t) {
        native.read = std::make_unique<detail::ReadRecorder>(p, id);
    }
};

} // namespace nustar::star2711b::solo
