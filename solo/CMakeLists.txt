add_library(nustar_solo_kernels STATIC
    solo/annotation/alignment.cpp
    solo/annotation/classify.cpp
    solo/annotation/exon_over_intron.cpp
    solo/annotation/exon_overlap.cpp
    solo/annotation/gene_overlap.cpp
    solo/annotation/index.cpp
    solo/barcode/extract.cpp
    solo/barcode/layout.cpp
    solo/barcode/match.cpp
    solo/barcode/match_options.cpp
    solo/barcode/read_barcode.cpp
    solo/barcode/sequence.cpp
    solo/barcode/whitelist.cpp
    solo/configuration/context.cpp
    solo/configuration/options.cpp
    solo/counting/counter.cpp
    solo/counting/feature_counter.cpp
    solo/counting/feature_recorder.cpp
    solo/counting/gene_umi.cpp
    solo/counting/input_record.cpp
    solo/counting/merge_workers.cpp
    solo/counting/multimapper_options.cpp
    solo/counting/process.cpp
    solo/counting/read_recorder.cpp
    solo/counting/read_records.cpp
    solo/counting/record_feature.cpp
    solo/counting/record_read.cpp
    solo/counting/redistribute.cpp
    solo/counting/smartseq.cpp
    solo/counting/transcript.cpp
    solo/counting/velocyto.cpp
    solo/filtering/cell_filter.cpp
    solo/filtering/empty_drops.cpp
    solo/filtering/options.cpp
    solo/matrix/load.cpp
    solo/matrix/statistics.cpp
    solo/matrix/write.cpp
    solo/support/binary_search.cpp
    solo/support/error.cpp
    solo/support/memory.cpp
    solo/support/streams.cpp
    solo/support/time.cpp
    solo/umi/deduplicate.cpp
    solo/umi/graph.cpp
    solo/umi/options.cpp
)
target_include_directories(nustar_solo_kernels PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_features(nustar_solo_kernels PRIVATE cxx_std_20)
set_target_properties(nustar_solo_kernels PROPERTIES CXX_EXTENSIONS OFF)
set_target_properties(nustar_solo_kernels PROPERTIES INTERPROCEDURAL_OPTIMIZATION_RELEASE OFF)
target_link_libraries(nustar_solo_kernels PUBLIC Threads::Threads)
if(OpenMP_CXX_FOUND)
    target_link_libraries(nustar_solo_kernels PUBLIC OpenMP::OpenMP_CXX)
endif()
add_library(nustar_solo STATIC solo/session.cpp solo/worker.cpp solo/filter.cpp)
target_link_libraries(nustar_solo PUBLIC nustar_contracts nustar_stats PRIVATE nustar_solo_kernels)
