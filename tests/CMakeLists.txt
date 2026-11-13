# Included from the project root; paths and binary names remain root-relative.
if(BUILD_TESTING)
    find_package(Python3 COMPONENTS Interpreter REQUIRED)
    add_executable(selection_contract_test tests/selection/contracts.cpp)
    target_link_libraries(selection_contract_test PRIVATE nustar_single_pipeline
        nustar_paired_pipeline nustar_defaults nustar_input)
    add_test(NAME selection_contract_test COMMAND selection_contract_test)
    add_executable(junction_strand_test tests/annotation/junction_strand.cpp)
    target_link_libraries(junction_strand_test PRIVATE nustar_annotation)
    add_test(NAME junction_strand_test COMMAND junction_strand_test)
    add_test(NAME star_cli COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/star_cli.py $<TARGET_FILE:nustar>)
    add_executable(single_projection_test tests/sam/single_projection.cpp)
    target_link_libraries(single_projection_test PRIVATE nustar_sam nustar_input)
    add_test(NAME single_projection_test COMMAND single_projection_test)
    add_executable(base_compare_test tests/seed_search/base_compare.cpp)
    target_link_libraries(base_compare_test PRIVATE nustar_contracts)
    add_test(NAME base_compare_test COMMAND base_compare_test)
    add_executable(sha256_correctness tests/common/sha256_correctness.cpp)
    target_link_libraries(sha256_correctness PRIVATE nustar_utils)
    add_test(NAME sha256_correctness COMMAND sha256_correctness)
    add_executable(sha256_portable_correctness tests/common/sha256_correctness.cpp common/sha256.cpp)
    target_link_libraries(sha256_portable_correctness PRIVATE nustar_contracts)
    target_compile_definitions(sha256_portable_correctness PRIVATE NUSTAR_SHA256_PORTABLE=1)
    add_test(NAME sha256_portable_correctness COMMAND sha256_portable_correctness)
    add_executable(window_workspace_reuse tests/window/workspace_reuse.cpp)
    target_link_libraries(window_workspace_reuse PRIVATE nustar_window)
    add_test(NAME window_workspace_reuse COMMAND window_workspace_reuse)
    add_executable(stage_progress_smoke tests/profiling/stage_progress_smoke.cpp)
    target_link_libraries(stage_progress_smoke PRIVATE nustar_timing)
    add_test(NAME stage_progress_smoke COMMAND stage_progress_smoke
        ${CMAKE_CURRENT_BINARY_DIR}/stage_progress_smoke.tsv)
    add_executable(packed_suffix_array tests/index/packed_suffix_array.cpp)
    target_link_libraries(packed_suffix_array PRIVATE nustar_sa_index)
    add_test(NAME packed_suffix_array COMMAND packed_suffix_array)
    add_executable(libsais40_correctness tests/index/libsais40.cpp)
    target_link_libraries(libsais40_correctness PRIVATE nustar_genomegenerate)
    add_test(NAME libsais40_correctness COMMAND libsais40_correctness)
    # Exercise large-alphabet 40-bit recursion without allocating >2^31 entries.
    add_executable(libsais40_compact_correctness tests/index/libsais40.cpp index/libsais40.cpp)
    target_compile_definitions(libsais40_compact_correctness PRIVATE NUSTAR_LIBSAIS40_FORCE_COMPACT_RECURSION=1)
    target_link_libraries(libsais40_compact_correctness PRIVATE nustar_genomegenerate)
    add_test(NAME libsais40_compact_correctness COMMAND libsais40_compact_correctness)
    add_executable(star_rank_correctness tests/index/star_rank.cpp)
    target_link_libraries(star_rank_correctness PRIVATE nustar_genomegenerate)
    add_test(NAME star_rank_correctness COMMAND star_rank_correctness)

    add_executable(genome_generate_smoke tests/index/genome_generate_smoke.cpp)
    target_link_libraries(genome_generate_smoke PRIVATE nustar_sa_index nustar_reference_service)
    add_test(NAME genome_generate_smoke COMMAND genome_generate_smoke)
    add_executable(seed_binding_smoke tests/seed_search/seed_binding_smoke.cpp)
    target_link_libraries(seed_binding_smoke PRIVATE nustar_seed nustar_reference_service Threads::Threads)
    add_test(NAME seed_binding_smoke COMMAND seed_binding_smoke
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/star_rank_probe.fa
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/star_rank_probe.gtf)
    add_executable(input_dispatch_smoke tests/reads/input_dispatch_smoke.cpp)
    target_link_libraries(input_dispatch_smoke PRIVATE nustar_input)
    add_test(NAME input_dispatch_smoke COMMAND input_dispatch_smoke)
    add_executable(fastq_io_smoke tests/reads/fastq_io_smoke.cpp)
    target_link_libraries(fastq_io_smoke PRIVATE nustar_input)
    add_test(NAME fastq_io_smoke COMMAND fastq_io_smoke
        ${CMAKE_CURRENT_BINARY_DIR}/fastq_io_smoke.data)

    add_executable(mapping_input_smoke tests/reads/mapping_input_smoke.cpp)
    target_link_libraries(mapping_input_smoke PRIVATE nustar_contracts)
    add_test(NAME mapping_input_smoke COMMAND mapping_input_smoke)

    add_executable(se_link_smoke tests/pipeline/se_link_smoke.cpp)
    target_link_libraries(se_link_smoke PRIVATE nustar_se)
    add_test(NAME se_link_smoke COMMAND se_link_smoke)

    add_executable(pe_link_smoke tests/pipeline/pe_link_smoke.cpp)
    target_link_libraries(pe_link_smoke PRIVATE nustar_pe)
    add_test(NAME pe_link_smoke COMMAND pe_link_smoke)

    add_test(NAME module_headers_and_boundaries
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/architecture/check_architecture.py
            ${CMAKE_CXX_COMPILER}
            ${CMAKE_CURRENT_SOURCE_DIR})

    add_test(NAME framework_help COMMAND nustar --help)

    set(NUSTAR_STAR_BASELINE "" CACHE FILEPATH
        "STAR 2.7.11b executable for the optional GenomeGenerate differential test")
    set(NUSTAR_STAR_DIFF_FASTA "${CMAKE_CURRENT_SOURCE_DIR}/../data/test/test.fa" CACHE FILEPATH
        "FASTA fixture for the optional GenomeGenerate differential test")
    set(NUSTAR_STAR_DIFF_GTF "${CMAKE_CURRENT_SOURCE_DIR}/../data/test/test.gtf" CACHE FILEPATH
        "GTF fixture for the optional GenomeGenerate differential test")
    set(NUSTAR_STAR_SEED_ORACLE "" CACHE FILEPATH
        "Instrumented STAR executable built by tools/oracle/build_star_seed_oracle.py")
    set(NUSTAR_STAR_STITCH_ORACLE "" CACHE FILEPATH
        "Instrumented STAR executable built by tools/oracle/build_star_stitch_oracle.py")
    set(NUSTAR_STAR_SEED_STAR_INDEX "" CACHE PATH
        "STAR index for the optional seed search differential test")
    set(NUSTAR_STAR_SEED_NUSTAR_INDEX "" CACHE PATH
        "NuSTAR index for the optional seed search differential test")
    if(NUSTAR_STAR_SEED_ORACLE AND NUSTAR_STAR_SEED_STAR_INDEX AND
       NUSTAR_STAR_SEED_NUSTAR_INDEX)
        add_test(NAME seed_search_star_diff
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/run_seed_diff.py
                --star-oracle ${NUSTAR_STAR_SEED_ORACLE}
                --nustar-trace $<TARGET_FILE:seed-trace>
                --star-index ${NUSTAR_STAR_SEED_STAR_INDEX}
                --nustar-index ${NUSTAR_STAR_SEED_NUSTAR_INDEX}
                --mate1 ${CMAKE_CURRENT_SOURCE_DIR}/../data/test/test_paired1.fastq
                --mate2 ${CMAKE_CURRENT_SOURCE_DIR}/../data/test/test_paired2.fastq
                --sample 1000)
        set_tests_properties(seed_search_star_diff PROPERTIES LABELS integration)
    endif()
    if(NUSTAR_STAR_BASELINE AND NUSTAR_STAR_SEED_STAR_INDEX AND
       NUSTAR_STAR_SEED_NUSTAR_INDEX)
        add_test(NAME single_end_sam_star_diff
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/benchmark/benchmark_se_sam.py
                --star ${NUSTAR_STAR_BASELINE}
                --nustar $<TARGET_FILE:se-sam-trace>
                --compare-sj
                --compare-stats
                --star-index ${NUSTAR_STAR_SEED_STAR_INDEX}
                --nustar-index ${NUSTAR_STAR_SEED_NUSTAR_INDEX}
                --reads ${CMAKE_CURRENT_SOURCE_DIR}/../docs/examples/STAR_2.7.11b_single_stitch_1000_2026-10-09/reads.fastq
                --work-dir ${CMAKE_CURRENT_BINARY_DIR}/single_end_sam_diff
                --repeats 1
                --output ${CMAKE_CURRENT_BINARY_DIR}/single_end_sam_diff.json)
        set_tests_properties(single_end_sam_star_diff PROPERTIES LABELS integration)
        add_test(NAME product_single_end_sam_star_diff
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/benchmark/benchmark_se_sam.py
                --star ${NUSTAR_STAR_BASELINE}
                --nustar $<TARGET_FILE:nustar>
                --nustar-mode product
                --compare-sj
                --compare-stats
                --star-index ${NUSTAR_STAR_SEED_STAR_INDEX}
                --nustar-index ${NUSTAR_STAR_SEED_NUSTAR_INDEX}
                --reads ${CMAKE_CURRENT_SOURCE_DIR}/../docs/examples/STAR_2.7.11b_single_stitch_1000_2026-10-09/reads.fastq
                --work-dir ${CMAKE_CURRENT_BINARY_DIR}/product_single_end_sam_diff
                --repeats 1
                --output ${CMAKE_CURRENT_BINARY_DIR}/product_single_end_sam_diff.json)
        set_tests_properties(product_single_end_sam_star_diff PROPERTIES LABELS integration)
        add_test(NAME product_single_end_unmapped_star_diff
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/benchmark/benchmark_se_sam.py
                --star ${NUSTAR_STAR_BASELINE}
                --nustar $<TARGET_FILE:nustar>
                --nustar-mode product
                --compare-sj
                --compare-stats
                --star-index ${NUSTAR_STAR_SEED_STAR_INDEX}
                --nustar-index ${NUSTAR_STAR_SEED_NUSTAR_INDEX}
                --reads ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/unmapped_read.fastq
                --work-dir ${CMAKE_CURRENT_BINARY_DIR}/product_single_end_unmapped_diff
                --repeats 1
                --output ${CMAKE_CURRENT_BINARY_DIR}/product_single_end_unmapped_diff.json)
        set_tests_properties(product_single_end_unmapped_star_diff PROPERTIES LABELS integration)
        add_test(NAME product_paired_sam_star_diff
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/benchmark/benchmark_se_sam.py
                --star ${NUSTAR_STAR_BASELINE}
                --nustar $<TARGET_FILE:nustar>
                --nustar-mode paired
                --compare-sj
                --compare-stats
                --star-index ${NUSTAR_STAR_SEED_STAR_INDEX}
                --nustar-index ${NUSTAR_STAR_SEED_NUSTAR_INDEX}
                --reads ${CMAKE_CURRENT_SOURCE_DIR}/../docs/examples/STAR_2.7.11b_paired_1000_2026-10-09/mate1.fastq
                --mate2 ${CMAKE_CURRENT_SOURCE_DIR}/../docs/examples/STAR_2.7.11b_paired_1000_2026-10-09/mate2.fastq
                --work-dir ${CMAKE_CURRENT_BINARY_DIR}/product_paired_sam_diff
                --repeats 1
                --output ${CMAKE_CURRENT_BINARY_DIR}/product_paired_sam_diff.json)
        set_tests_properties(product_paired_sam_star_diff PROPERTIES LABELS integration)
    endif()
    if(NUSTAR_STAR_STITCH_ORACLE AND NUSTAR_STAR_SEED_STAR_INDEX AND
       NUSTAR_STAR_SEED_NUSTAR_INDEX)
        add_test(NAME paired_window_stitch_star_diff
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/run_pe_window_diff.py
                --star-oracle ${NUSTAR_STAR_STITCH_ORACLE}
                --nustar-trace $<TARGET_FILE:window-trace>
                --stitch-trace $<TARGET_FILE:paired-stitch-trace>
                --star-index ${NUSTAR_STAR_SEED_STAR_INDEX}
                --nustar-index ${NUSTAR_STAR_SEED_NUSTAR_INDEX}
                --mate1 ${CMAKE_CURRENT_SOURCE_DIR}/../data/test/test_paired1.fastq
                --mate2 ${CMAKE_CURRENT_SOURCE_DIR}/../data/test/test_paired2.fastq
                --sample 1000)
        set_tests_properties(paired_window_stitch_star_diff PROPERTIES LABELS integration)
        add_test(NAME single_stitch_star_diff
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/run_stitch_diff.py
                --star-oracle ${NUSTAR_STAR_STITCH_ORACLE}
                --nustar-trace $<TARGET_FILE:stitch-trace>
                --star-index ${NUSTAR_STAR_SEED_STAR_INDEX}
                --nustar-index ${NUSTAR_STAR_SEED_NUSTAR_INDEX}
                --reads ${CMAKE_CURRENT_SOURCE_DIR}/../data/test/test_paired1.fastq
                --sample 1000)
        set_tests_properties(single_stitch_star_diff PROPERTIES LABELS integration)
    endif()
    if(NUSTAR_STAR_STITCH_ORACLE AND NUSTAR_STAR_BASELINE)
        add_test(NAME single_stitch_star_edge_diff
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/run_stitch_edge_diff.py
                --star ${NUSTAR_STAR_BASELINE}
                --nustar $<TARGET_FILE:nustar>
                --compare-index $<TARGET_FILE:compare-genome-indices>
                --star-oracle ${NUSTAR_STAR_STITCH_ORACLE}
                --nustar-trace $<TARGET_FILE:stitch-trace>
                --fasta ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/star_rank_probe.fa
                --gtf ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/star_rank_probe.gtf
                --reads ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/seed_edge_reads.fastq)
        set_tests_properties(single_stitch_star_edge_diff PROPERTIES LABELS integration)
    endif()
    if(NUSTAR_STAR_SEED_ORACLE AND NUSTAR_STAR_BASELINE)
        add_test(NAME seed_search_star_edge_diff
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/run_seed_edge_diff.py
                --star ${NUSTAR_STAR_BASELINE}
                --nustar $<TARGET_FILE:nustar>
                --star-oracle ${NUSTAR_STAR_SEED_ORACLE}
                --nustar-trace $<TARGET_FILE:seed-trace>
                --fasta ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/star_rank_probe.fa
                --gtf ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/star_rank_probe.gtf
                --reads ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/seed_edge_reads.fastq)
        set_tests_properties(seed_search_star_edge_diff PROPERTIES LABELS integration)
    endif()
    if(NUSTAR_STAR_BASELINE)
        add_test(NAME genome_generate_star_diff
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/run_genome_generate_diff.py
                --star ${NUSTAR_STAR_BASELINE}
                --nustar $<TARGET_FILE:nustar>
                --compare-index $<TARGET_FILE:compare-genome-indices>
                --fasta ${NUSTAR_STAR_DIFF_FASTA}
                --gtf ${NUSTAR_STAR_DIFF_GTF}
                --overhang 150 --chr-bin-bits 18 --prefix-bases 10)
        set_tests_properties(genome_generate_star_diff PROPERTIES LABELS integration)
        add_test(NAME genome_generate_star_diff_no_gtf
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/run_genome_generate_diff.py
                --star ${NUSTAR_STAR_BASELINE}
                --nustar $<TARGET_FILE:nustar>
                --compare-index $<TARGET_FILE:compare-genome-indices>
                --fasta ${NUSTAR_STAR_DIFF_FASTA}
                --chr-bin-bits 18 --prefix-bases 10)
        set_tests_properties(genome_generate_star_diff_no_gtf PROPERTIES LABELS integration)
        add_test(NAME genome_generate_star_rank_boundaries
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/run_genome_generate_diff.py
                --star ${NUSTAR_STAR_BASELINE}
                --nustar $<TARGET_FILE:nustar>
                --compare-index $<TARGET_FILE:compare-genome-indices>
                --fasta ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/star_rank_probe.fa
                --gtf ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/star_rank_probe.gtf
                --overhang 16 --chr-bin-bits 8 --prefix-bases 4)
        set_tests_properties(genome_generate_star_rank_boundaries PROPERTIES LABELS integration)
    endif()
endif()

# Native 2-pass/Solo compatibility suites construct their own small indices.
if(BUILD_TESTING)
    set(NUSTAR_FEATURE_STAR_BASELINE "" CACHE FILEPATH
        "STAR 2.7.11b executable for native 2-pass and Solo differential suites")
    if(NUSTAR_FEATURE_STAR_BASELINE)
        add_test(NAME junction_strand_star_diff COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/junction_strand.py
            $<TARGET_FILE:nustar> ${NUSTAR_FEATURE_STAR_BASELINE}
            ${CMAKE_CURRENT_BINARY_DIR}/junction-strand-diff)
        set_tests_properties(junction_strand_star_diff PROPERTIES LABELS integration)
        add_test(NAME paired_overlap_star_diff COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/paired_overlap.py
            $<TARGET_FILE:nustar> ${NUSTAR_FEATURE_STAR_BASELINE}
            ${CMAKE_CURRENT_BINARY_DIR}/paired-overlap-diff)
        set_tests_properties(paired_overlap_star_diff PROPERTIES LABELS integration)
        add_test(NAME junction_padding_star_diff COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/junction_padding.py
            $<TARGET_FILE:nustar> ${NUSTAR_FEATURE_STAR_BASELINE}
            $<TARGET_FILE:compare-genome-indices>
            ${CMAKE_CURRENT_BINARY_DIR}/junction-padding-diff)
        set_tests_properties(junction_padding_star_diff PROPERTIES LABELS integration)
        add_test(NAME two_pass_star_diff COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/two_pass.py
            $<TARGET_FILE:nustar> ${NUSTAR_FEATURE_STAR_BASELINE}
            ${CMAKE_CURRENT_BINARY_DIR}/two-pass-diff)
        add_test(NAME solo_star_diff COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/solo.py
            $<TARGET_FILE:nustar> ${NUSTAR_FEATURE_STAR_BASELINE}
            ${CMAKE_CURRENT_BINARY_DIR}/solo-diff)
        add_test(NAME solo_modes_star_diff COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/solo_modes.py
            $<TARGET_FILE:nustar> ${NUSTAR_FEATURE_STAR_BASELINE}
            ${CMAKE_CURRENT_BINARY_DIR}/solo-diff)
        set_tests_properties(solo_star_diff PROPERTIES FIXTURES_SETUP solo_reference)
        set_tests_properties(solo_modes_star_diff PROPERTIES FIXTURES_REQUIRED solo_reference)
        add_test(NAME bam_star_diff COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/bam.py
            $<TARGET_FILE:nustar> ${NUSTAR_FEATURE_STAR_BASELINE}
            ${CMAKE_CURRENT_BINARY_DIR}/solo-diff ${CMAKE_CURRENT_BINARY_DIR}/bam-diff)
        set_tests_properties(bam_star_diff PROPERTIES FIXTURES_REQUIRED solo_reference
            FIXTURES_SETUP bam_reference LABELS integration)
        add_test(NAME bam_pair_flags_star_diff COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/bam_pair_flags.py
            $<TARGET_FILE:nustar> ${NUSTAR_FEATURE_STAR_BASELINE}
            ${CMAKE_CURRENT_BINARY_DIR}/solo-diff ${CMAKE_CURRENT_BINARY_DIR}/bam-diff)
        set_tests_properties(bam_pair_flags_star_diff PROPERTIES
            FIXTURES_REQUIRED "solo_reference;bam_reference" LABELS integration)
        set_tests_properties(two_pass_star_diff solo_star_diff solo_modes_star_diff
            PROPERTIES LABELS integration)
    endif()
endif()
