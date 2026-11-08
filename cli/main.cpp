#include "genome_generate.hpp"
#include "align_reads.hpp"
#include "solo_filter.hpp"
#include <iostream>
#include <string_view>

int main(int argc, char** argv) {
    if (argc == 2 && std::string_view(argv[1]) == "--help") {
        std::cout << "NuSTAR\n"
                     "Usage: nustar --runMode genomeGenerate --genomeFastaFiles FILE... --genomeDir DIR [--sjdbGTFfile FILE --sjdbOverhang N] [--genomeChrBinNbits N] [--genomeSAindexNbases N] [--genomeSAalgorithm libsais64|libsais40] [--runThreadN N] [--profileGenomeGenerate]\n"
                     "Upgrade an existing index without rebuilding SA: upgradeIndex SOURCE_DIR NEW_DIR.\n"
                     "Align reads: [--runMode alignReads] --genomeDir DIR --readFilesIn MATE1.fastq [MATE2.fastq] --runThreadN N --outFileNamePrefix PREFIX\n"
                     "BAM input: --readFilesType SAM SE|PE; --readFilesSAMattrKeep All|None|tags.\n"
                     "BAM output: --outBAMcompression -1..9; --limitBAMsortRAM bytes; --outSAMattributes tags.\n"
                     "Output: --outSAMtype SAM|None or BAM Unsorted|SortedByCoordinate (both allowed); --outSAMunmapped None|Within (default None); --outSAMorder Paired|PairedKeepInputOrder.\n"
                     "Two-pass: --twopassMode Basic [--twopass1readsN N] [--sjdbInsertSave Basic|All].\n"
                     "Overlapping mates: --peOverlapNbasesMin N (default 0, disabled); --peOverlapMMp fraction (default 0.01).\n"
                     "Solo: --soloType CB_UMI_Simple|CB_UMI_Complex|SmartSeq; supports the --solo* parameters listed in readme.md.\n"
                     "Place the barcode FASTQ last in --readFilesIn, or use --soloBarcodeMate with fixed clipping.\n"
                     "SmartSeq: --readFilesManifest FILE --soloType SmartSeq --soloUMIdedup Exact.\n"
                     "Filter an existing matrix: --runMode soloCellFiltering RAW_DIR OUTPUT_PREFIX [--soloCellFilter ...].\n"
                     "NuSTAR options: --io-backend mmap|buffered, --profile-stages (single-end), --profile-runtime.\n"
                     "See readme.md for examples and supported options; CODE_STRUCTURE.md for the source layout.\n"
                     "The generated reference.bin and index.sa use NuSTAR reference v1 / packed SA index v4 formats.\n";
        return 0;
    }
    if (argc>=2 && std::string_view(argv[1])=="upgradeIndex")
        return runUpgradeIndex(argc,argv);
    if (argc == 2 && std::string_view(argv[1]) == "--version") {
        std::cout << "NuSTAR 0.2.0\n";
        return 0;
    }
    for(int i=1;i+1<argc;++i)
        if(std::string_view(argv[i])=="--runMode" && std::string_view(argv[i+1])=="soloCellFiltering")
            return runSoloCellFiltering(argc,argv);
    for (int i = 1; i + 1 < argc; ++i)
        if (std::string_view(argv[i]) == "--runMode" &&
            std::string_view(argv[i+1]) == "genomeGenerate")
            return runGenomeGenerate(argc, argv, 1);
    if (argc >= 2 && std::string_view(argv[1]).starts_with("--"))
        return runAlignReads(argc, argv);
    std::cerr << "NuSTAR: unsupported command; see --help.\n";
    return 2;
}
