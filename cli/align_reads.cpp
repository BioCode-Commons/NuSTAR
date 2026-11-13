#include "align_reads.hpp"
#include "../runtime/paired_alignment.hpp"
#include "../runtime/single_alignment.hpp"
#include "../sam/config.hpp"
#include "align_paired.hpp"
#include "align_single.hpp"

#include <charconv>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <vector>

int runAlignReads(int argc, char** argv) {
    try {
        std::string genomeDir = "./GenomeDir/", outFileNamePrefix = "./";
        std::vector<std::string> readFilesIn;
        std::optional<std::filesystem::path> manifest;
        std::uint32_t runThreadN = 1;
        std::string backend = "pread";
        bool within = false, profile = false, runtime_profile = false, bam_paired = false;
        nustar::star2711b::RunControl control;
        auto solo_config = std::make_shared<nustar::star2711b::solo::Config>();
        for (int i = 1; i < argc;) {
            const std::string_view flag(argv[i++]);
            auto value = [&]() -> std::string_view {
                if (i == argc || std::string_view(argv[i]).starts_with("--"))
                    throw std::invalid_argument("missing value for " + std::string(flag));
                return argv[i++];
            };
            if (flag == "--runMode") {
                if (value() != "alignReads")
                    throw std::invalid_argument("unsupported --runMode");
            } else if (flag == "--genomeDir")
                genomeDir = value();
            else if (flag == "--readFilesIn") {
                readFilesIn.clear();
                while (i < argc && !std::string_view(argv[i]).starts_with("--"))
                    readFilesIn.emplace_back(argv[i++]);
                if (readFilesIn.empty() || readFilesIn.size() > 3)
                    throw std::invalid_argument("--readFilesIn requires RNA FASTQ paths and an "
                                                "optional Solo barcode FASTQ");
                for (const auto& path : readFilesIn)
                    if (path.find(',') != std::string::npos)
                        throw std::invalid_argument(
                            "comma-separated read file lists are not implemented");
            } else if (flag == "--readFilesManifest")
                manifest = std::string(value());
            else if (flag == "--clip5pNbases" || flag == "--clip3pNbases") {
                const unsigned end = flag == "--clip5pNbases" ? 0 : 1;
                unsigned mate = 0;
                do {
                    const auto token = value();
                    std::uint64_t amount;
                    const auto parse =
                        std::from_chars(token.data(), token.data() + token.size(), amount);
                    if (parse.ec != std::errc{} || parse.ptr != token.data() + token.size() ||
                        mate >= 2)
                        throw std::invalid_argument("invalid clipping value");
                    control.clipping[mate++][end] = amount;
                } while (i < argc && !std::string_view(argv[i]).starts_with("--"));
                if (mate == 1 && control.clipping[0][end] == 0)
                    control.clipping[1][end] = 0;
            } else if (flag == "--outFileNamePrefix")
                outFileNamePrefix = value();
            else if (flag == "--runThreadN") {
                const auto token = value();
                const auto parsed =
                    std::from_chars(token.data(), token.data() + token.size(), runThreadN);
                if (parsed.ec != std::errc{} || parsed.ptr != token.data() + token.size() ||
                    runThreadN == 0 || runThreadN > 256)
                    throw std::invalid_argument("--runThreadN must be between 1 and 256");
            } else if (flag == "--peOverlapNbasesMin") {
                const auto token=value();
                const auto parsed=std::from_chars(token.data(),token.data()+token.size(),control.overlap.peOverlapNbasesMin);
                if (parsed.ec!=std::errc{} || parsed.ptr!=token.data()+token.size())
                    throw std::invalid_argument("invalid --peOverlapNbasesMin");
            } else if (flag == "--peOverlapMMp") {
                const auto token=value();
                const auto parsed=std::from_chars(token.data(),token.data()+token.size(),control.overlap.peOverlapMMp);
                if (parsed.ec!=std::errc{} || parsed.ptr!=token.data()+token.size() ||
                    !std::isfinite(control.overlap.peOverlapMMp) || control.overlap.peOverlapMMp<0 || control.overlap.peOverlapMMp>=1)
                    throw std::invalid_argument("--peOverlapMMp must be >=0 and <1");
            } else if (flag == "--outSAMtype") {
                const auto type = value();
                if (type != "SAM" && type != "None" && type != "BAM")
                    throw std::invalid_argument("--outSAMtype supports SAM, BAM, or None");
                control.sam_enabled = type != "None";
                control.bam.unsorted = control.bam.coordinate = false;
                if (type == "BAM") {
                    do {
                        const auto order = value();
                        if (order == "Unsorted")
                            control.bam.unsorted = true;
                        else if (order == "SortedByCoordinate")
                            control.bam.coordinate = true;
                        else
                            throw std::invalid_argument(
                                "BAM output requires Unsorted and/or SortedByCoordinate");
                    } while (i < argc && !std::string_view(argv[i]).starts_with("--"));
                }
            } else if (flag == "--readFilesType") {
                const auto type = value();
                if (type == "SAM") {
                    control.bam_input = true;
                    const auto layout = value();
                    if (layout != "SE" && layout != "PE")
                        throw std::invalid_argument("--readFilesType SAM requires SE or PE");
                    bam_paired = layout == "PE";
                } else if (type != "Fastx")
                    throw std::invalid_argument("--readFilesType supports Fastx or SAM SE/PE (SAM "
                                                "or BAM detected automatically)");
            } else if (flag == "--readFilesSAMattrKeep") {
                control.input_attributes.clear();
                do {
                    control.input_attributes.emplace_back(value());
                } while (i < argc && !std::string_view(argv[i]).starts_with("--"));
            } else if (flag == "--outSAMattributes") {
                control.output_attributes.clear();
                do {
                    const std::string tag(value());
                    if (tag == "Standard")
                        for (const auto* name : {"NH", "HI", "AS", "nM"})
                            control.output_attributes.emplace_back(name);
                    else
                        control.output_attributes.push_back(tag);
                } while (i < argc && !std::string_view(argv[i]).starts_with("--"));
            } else if (flag == "--outBAMcompression") {
                control.bam.compression = std::stoi(std::string(value()));
                if (control.bam.compression < -1 || control.bam.compression > 9)
                    throw std::invalid_argument("--outBAMcompression must be -1 through 9");
            } else if (flag == "--limitBAMsortRAM") {
                const auto token = value();
                const auto parsed = std::from_chars(token.data(), token.data() + token.size(),
                                                    control.bam.sort_memory);
                if (parsed.ec != std::errc{} || parsed.ptr != token.data() + token.size())
                    throw std::invalid_argument("invalid --limitBAMsortRAM");
                if (control.bam.sort_memory == 0)
                    control.bam.sort_memory = 256ULL << 20;
            } else if (flag == "--outSAMunmapped") {
                const auto mode = value();
                if (mode != "None" && mode != "Within")
                    throw std::invalid_argument("--outSAMunmapped supports None or Within");
                within = mode == "Within";
            } else if (flag == "--outSAMorder") {
                if (value() != "Paired")
                    throw std::invalid_argument(
                        "--outSAMorder only supports Paired; input-order output has been removed");
            } else if (flag == "--io-backend") {
                backend = value();
                if (backend != "pread" && backend != "mmap" && backend != "buffered")
                    throw std::invalid_argument("--io-backend supports pread (mmap and buffered are compatibility aliases)");
            } else if (flag == "--twopassMode") {
                const auto mode = value();
                if (mode != "None" && mode != "Basic")
                    throw std::invalid_argument("--twopassMode supports None or Basic");
                control.twopass = mode == "Basic";
            } else if (flag == "--sjdbInsertSave") {
                const auto mode = value();
                if (mode != "Basic" && mode != "All")
                    throw std::invalid_argument("--sjdbInsertSave supports Basic or All");
                control.save_inserted_index = mode == "All";
            } else if (flag == "--twopass1readsN" || flag == "--readMapNumber" ||
                       flag == "--sjdbOverhang" || flag == "--limitSjdbInsertNsj") {
                const auto token = value();
                std::uint64_t number = 0;
                if (token == "-1" && (flag == "--twopass1readsN" || flag == "--readMapNumber"))
                    number = UINT64_MAX;
                else {
                    const auto parsed =
                        std::from_chars(token.data(), token.data() + token.size(), number);
                    if (parsed.ec != std::errc{} || parsed.ptr != token.data() + token.size())
                        throw std::invalid_argument("invalid " + std::string(flag));
                }
                if (flag == "--twopass1readsN")
                    control.pass1_reads = number;
                else if (flag == "--readMapNumber")
                    control.read_limit = number;
                else if (flag == "--sjdbOverhang")
                    control.sjdb_overhang = number;
                else
                    control.junction_limit = number;
            } else if (flag == "--profile-stages")
                profile = true;
            else if (flag == "--profile-runtime")
                runtime_profile = true;
            else if (flag.starts_with("--solo")) {
                solo_config->arguments.emplace_back(flag);
                auto first = value();
                solo_config->arguments.emplace_back(first);
                if (flag == "--soloType")
                    solo_config->type = first;
                if (flag == "--soloBarcodeMate") {
                    const auto parsed = std::from_chars(first.data(), first.data() + first.size(),
                                                        solo_config->barcode_mate);
                    if (parsed.ec != std::errc{} || parsed.ptr != first.data() + first.size())
                        throw std::invalid_argument("invalid --soloBarcodeMate");
                }
                while (i < argc && !std::string_view(argv[i]).starts_with("--"))
                    solo_config->arguments.emplace_back(argv[i++]);
            } else
                throw std::invalid_argument("unsupported alignment option: " + std::string(flag));
        }
        control.bam.threads = runThreadN;
        solo_config->clipping = control.clipping;
        solo_config->bam_input = control.bam_input;
        solo_config->coordinate_bam = control.bam.coordinate;
        solo_config->output_attributes = control.output_attributes;
        for (const auto& tag : control.output_attributes) {
            if (nustar::star2711b::isSoloAttribute(tag) && solo_config->type == "None")
                throw std::invalid_argument("Solo output attributes require --soloType");
        }
        if (manifest) {
            if (solo_config->type != "SmartSeq")
                throw std::invalid_argument(
                    "--readFilesManifest currently requires --soloType SmartSeq");
            std::ifstream file(*manifest);
            if (!file)
                throw std::invalid_argument("cannot open read manifest");
            std::string line;
            bool paired = false;
            while (std::getline(file, line)) {
                if (line.find_first_not_of(" \t\r") == std::string::npos)
                    continue;
                const auto tab1 = line.find('\t'), tab2 = line.find('\t', tab1 + 1);
                if (tab1 == std::string::npos || tab2 == std::string::npos)
                    throw std::invalid_argument("manifest requires three tab-separated columns");
                const auto mate1 = line.substr(0, tab1),
                           mate2 = line.substr(tab1 + 1, tab2 - tab1 - 1);
                const bool row_paired = mate2 != "-";
                if (!control.input_groups.empty() && row_paired != paired)
                    throw std::invalid_argument("manifest mixes SE and PE input");
                paired = row_paired;
                control.input_groups.push_back(
                    nustar::star2711b::InputGroup{{mate1, paired ? mate2 : ""}});
                auto group = line.substr(tab2 + 1);
                if (!group.empty() && group.back() == '\r')
                    group.pop_back();
                solo_config->read_group_lines.push_back(group.starts_with("ID:") ? group
                                                                                 : "ID:" + group);
                if (group.starts_with("ID:"))
                    group.erase(0, 3);
                group.resize(group.find_first_of("\t\r") == std::string::npos
                                 ? group.size()
                                 : group.find_first_of("\t\r"));
                solo_config->read_groups.push_back(group);
            }
            if (control.input_groups.empty())
                throw std::invalid_argument("empty read manifest");
            readFilesIn = {control.input_groups[0].mate_paths[0].string()};
            if (paired)
                readFilesIn.push_back(control.input_groups[0].mate_paths[1].string());
        }
        if (solo_config->type != "None") {
            control.solo_config = solo_config;
            if (!control.bam_input && solo_config->type != "SmartSeq" &&
                solo_config->barcode_mate == 0) {
                if (readFilesIn.size() < 2)
                    throw std::invalid_argument("Solo needs RNA and barcode reads");
                control.barcode_fastq = readFilesIn.back();
                readFilesIn.pop_back();
            }
        }
        if (control.bam_input) {
            if (readFilesIn.size() != 1 || manifest)
                throw std::invalid_argument("SAM/BAM input requires one file");
            if (bam_paired)
                readFilesIn.push_back(readFilesIn[0]);
        }
        if (solo_config->barcode_mate > readFilesIn.size())
            throw std::invalid_argument("--soloBarcodeMate exceeds RNA mate count");
        if (readFilesIn.size() > 2)
            throw std::invalid_argument("at most two RNA mates are supported");
        if (control.twopass && control.pass1_reads == 0)
            throw std::invalid_argument("--twopass1readsN must be >0 or -1");
        if (readFilesIn.empty())
            throw std::invalid_argument("--readFilesIn is required");
        if (profile && readFilesIn.size() == 2)
            throw std::invalid_argument(
                "--profile-stages currently supports single-end reads only");
        std::string command;
        for (int i = 0; i < argc; ++i) {
            if (i)
                command += ' ';
            command += argv[i];
        }
        const std::filesystem::path sam(outFileNamePrefix + "Aligned.out.sam");
        if (!sam.parent_path().empty())
            std::filesystem::create_directories(sam.parent_path());
        using namespace nustar::star2711b;
        if (readFilesIn.size() == 1)
            return runAlignSingleSam(SingleSamOptions{genomeDir, readFilesIn[0], sam,
                                                      std::move(command), profile, runThreadN,
                                                      backend, outFileNamePrefix, within,
                                                      control},
                                     runtime_profile);
        return runAlignPairedSam(PairedSamOptions{genomeDir, readFilesIn[0], readFilesIn[1], sam,
                                                  std::move(command), runThreadN, backend,
                                                  outFileNamePrefix, within, control},
                                 runtime_profile);
    } catch (const std::exception& error) {
        std::cerr << "alignReads: " << error.what() << '\n';
        return 2;
    }
}
