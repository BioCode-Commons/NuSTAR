find_package(ZLIB REQUIRED)
add_library(nustar_hts STATIC
    kfunc.c
    knetfile.c
    kstring.c
    bgzf.c
    faidx.c
    hfile.c
    hfile_net.c
    hts.c
    sam.c
    synced_bcf_reader.c
    vcf_sweep.c
    tbx.c
    vcf.c
    vcfutils.c
    cram/cram_codecs.c
    cram/cram_decode.c
    cram/cram_encode.c
    cram/cram_index.c
    cram/cram_io.c
    cram/cram_samtools.c
    cram/cram_stats.c
    cram/files.c
    cram/mFILE.c
    cram/md5.c
    cram/open_trace_file.c
    cram/pooled_alloc.c
    cram/sam_header.c
    cram/string_alloc.c
    cram/thread_pool.c
    cram/vlen.c
    cram/zfio.c
)
target_include_directories(nustar_hts PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(nustar_hts PRIVATE SAMTOOLS=1)
target_link_libraries(nustar_hts PUBLIC ZLIB::ZLIB Threads::Threads m)
set_target_properties(nustar_hts PROPERTIES C_STANDARD 99)
