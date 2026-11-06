cmake_minimum_required(VERSION 3.20)
project(NuSTAR VERSION 0.2.0 LANGUAGES C CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# The recursive stitching hot path crosses translation-unit boundaries.
# Release IPO lets the compiler optimize those calls while retaining a
# portable non-IPO build when the selected toolchain does not support it.
option(NUSTAR_ENABLE_RELEASE_IPO "Enable interprocedural optimization in Release" ON)
if(NUSTAR_ENABLE_RELEASE_IPO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT NUSTAR_IPO_SUPPORTED OUTPUT NUSTAR_IPO_ERROR
        LANGUAGES CXX)
    if(NUSTAR_IPO_SUPPORTED)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
    else()
        message(STATUS "Release IPO unavailable: ${NUSTAR_IPO_ERROR}")
    endif()
endif()

include(CTest)
find_package(Threads REQUIRED)
find_package(OpenMP QUIET COMPONENTS C CXX)

# Third-party implementation is isolated from the index adapters.
add_library(nustar_libsais STATIC third_party/libsais/libsais.c third_party/libsais/libsais64.c)
set_target_properties(nustar_libsais PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED ON)
target_include_directories(nustar_libsais PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/third_party/libsais)
if(OpenMP_C_FOUND)
    target_compile_definitions(nustar_libsais PUBLIC LIBSAIS_OPENMP)
    target_link_libraries(nustar_libsais PUBLIC OpenMP::OpenMP_C)
endif()

add_library(nustar_contracts INTERFACE)
target_compile_features(nustar_contracts INTERFACE cxx_std_20)
target_include_directories(nustar_contracts INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

add_library(nustar_utils STATIC common/sha256.cpp common/file_identity.cpp)
target_link_libraries(nustar_utils PUBLIC nustar_contracts)
add_library(nustar_timing STATIC profiling/stage_progress.cpp)
target_link_libraries(nustar_timing PUBLIC nustar_contracts)

add_library(nustar_reference STATIC reference/genome.cpp reference/synthetic_junction.cpp)
target_link_libraries(nustar_reference PUBLIC nustar_utils Threads::Threads)
add_library(nustar_annotation STATIC annotation/gtf_loader.cpp
    annotation/table_writer.cpp annotation/junction_database.cpp)
target_link_libraries(nustar_annotation PUBLIC nustar_contracts)
add_library(nustar_reference_service STATIC runtime/reference_service.cpp)
target_link_libraries(nustar_reference_service PUBLIC nustar_reference nustar_annotation)

add_library(nustar_genomegenerate STATIC index/sais.cpp
    index/wide_suffix_array.cpp index/star_rank.cpp)
target_link_libraries(nustar_genomegenerate PUBLIC nustar_contracts nustar_libsais)
if(OpenMP_CXX_FOUND)
    target_link_libraries(nustar_genomegenerate PUBLIC OpenMP::OpenMP_CXX)
endif()
add_library(nustar_sa_index STATIC index/suffix_array_index.cpp
    index/packed_suffix_array.cpp index/star_prefix.cpp index/physical_sequence.cpp index/extension.cpp index/junction_insert.cpp)
target_link_libraries(nustar_sa_index PUBLIC nustar_reference nustar_genomegenerate)

if(OpenMP_CXX_FOUND)
    target_link_libraries(nustar_sa_index PUBLIC OpenMP::OpenMP_CXX)
endif()

add_library(nustar_seed STATIC seed_search/star_search.cpp)
target_link_libraries(nustar_seed PUBLIC nustar_sa_index)
add_library(nustar_window STATIC window/star_window.cpp)
target_link_libraries(nustar_window PUBLIC nustar_contracts)
add_library(nustar_stitching STATIC stitching/single_stitch.cpp stitching/paired_stitch.cpp
    stitching/extend_alignment.cpp stitching/blocks_overlap.cpp
    stitching/stitch_segments.cpp stitching/stitch_window.cpp)
target_link_libraries(nustar_stitching PUBLIC nustar_reference)
add_library(nustar_selection STATIC selection/star_select.cpp)
target_link_libraries(nustar_selection PUBLIC nustar_contracts)

add_library(nustar_input STATIC reads/dispatch.cpp reads/fastq.cpp reads/fastq_prepare.cpp)
target_link_libraries(nustar_input PUBLIC nustar_contracts)
add_library(nustar_sam STATIC sam/encode.cpp sam/star_single.cpp sam/star_paired.cpp)
target_link_libraries(nustar_sam PUBLIC nustar_contracts)
add_library(nustar_sj STATIC splice_junction/star_junction.cpp splice_junction/table_writer.cpp)
target_link_libraries(nustar_sj PUBLIC nustar_contracts)
add_library(nustar_stats STATIC statistics/star_statistics.cpp)
target_link_libraries(nustar_stats PUBLIC nustar_contracts)
add_subdirectory(third_party/htslib)
add_library(nustar_bam STATIC bam/codec.cpp bam/source.cpp bam/sink.cpp)
target_link_libraries(nustar_bam PUBLIC nustar_output nustar_hts)

add_library(nustar_output STATIC output/file_sink.cpp)
target_link_libraries(nustar_output PUBLIC nustar_sam nustar_stats)

add_library(nustar_single_mapping STATIC pipeline/single_end/mapping_composition.cpp)
target_link_libraries(nustar_single_mapping PUBLIC nustar_seed nustar_window nustar_stitching)
add_library(nustar_paired_mapping STATIC pipeline/paired_end/mapping_composition.cpp)
target_link_libraries(nustar_paired_mapping PUBLIC nustar_seed nustar_window nustar_stitching)
add_library(nustar_single_pipeline STATIC pipeline/single_end/pipeline.cpp)
target_link_libraries(nustar_single_pipeline PUBLIC nustar_single_mapping
    nustar_selection nustar_sam nustar_sj nustar_stats)
add_library(nustar_paired_pipeline STATIC pipeline/paired_end/pipeline.cpp
    pipeline/paired_end/star_overlap.cpp)
target_link_libraries(nustar_paired_pipeline PUBLIC nustar_paired_mapping nustar_single_mapping
    nustar_selection nustar_sam nustar_sj nustar_stats)

add_library(nustar_parallel STATIC runtime/ordered_alignment.cpp)
target_link_libraries(nustar_parallel PUBLIC nustar_input nustar_output
    nustar_sj nustar_stats Threads::Threads)
add_library(nustar_defaults STATIC runtime/defaults.cpp)
target_link_libraries(nustar_defaults PUBLIC nustar_contracts)
# Keep the established product entry libraries and command-line interfaces.
add_library(nustar_se STATIC runtime/single_runner.cpp runtime/single_alignment.cpp)
target_link_libraries(nustar_se PUBLIC nustar_single_pipeline nustar_defaults nustar_parallel)
add_library(nustar_pe STATIC runtime/paired_runner.cpp runtime/paired_alignment.cpp)
target_link_libraries(nustar_pe PUBLIC nustar_paired_pipeline nustar_defaults nustar_parallel)
add_library(nustar_twopass STATIC runtime/two_pass.cpp)
target_link_libraries(nustar_twopass PUBLIC nustar_se nustar_pe nustar_reference_service nustar_sa_index)
add_library(nustar_runtime STATIC runtime/planner.cpp)
target_link_libraries(nustar_runtime PUBLIC nustar_se nustar_pe nustar_reference_service)

add_executable(nustar cli/main.cpp cli/genome_generate.cpp
    cli/align_single.cpp cli/align_paired.cpp cli/align_reads.cpp cli/solo_filter.cpp)
target_link_libraries(nustar PRIVATE nustar_runtime nustar_twopass nustar_sa_index)

include(tools/CMakeLists.txt)
include(tests/CMakeLists.txt)

include(solo/CMakeLists.txt)
target_link_libraries(nustar_se PUBLIC nustar_solo nustar_bam)
target_link_libraries(nustar_pe PUBLIC nustar_solo nustar_bam)
