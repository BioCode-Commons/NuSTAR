# Diagnostic executables remain available with BUILD_TESTING=OFF.
add_executable(compare-genome-indices tests/integration/compare_index.cpp)
target_link_libraries(compare-genome-indices PRIVATE nustar_sa_index)

add_executable(seed-trace tools/trace/seed_trace.cpp)
target_link_libraries(seed-trace PRIVATE nustar_seed nustar_timing)

add_executable(seed-benchmark tools/benchmark/seed_benchmark.cpp)
target_link_libraries(seed-benchmark PRIVATE nustar_seed nustar_utils)

add_executable(stitch-trace tools/trace/stitch_trace.cpp)
target_link_libraries(stitch-trace PRIVATE nustar_single_mapping nustar_timing)

add_executable(paired-stitch-trace tools/trace/paired_stitch_trace.cpp)
target_link_libraries(paired-stitch-trace PRIVATE nustar_paired_mapping nustar_timing)

add_executable(window-trace tools/trace/window_trace.cpp)
target_link_libraries(window-trace PRIVATE nustar_window nustar_seed)

add_executable(se-sam-trace tools/trace/se_sam_trace.cpp)
target_link_libraries(se-sam-trace PRIVATE nustar_se)

# Used only to supply textual SAM to the upstream STAR benchmark.
add_executable(bam-to-sam tools/benchmark/bam_to_sam.cpp)
target_link_libraries(bam-to-sam PRIVATE nustar_hts)
