#include "session.hpp"
#include "barcode/sequence.hpp"
#include "counting/counter.hpp"
#include "worker_state.hpp"
#include <fstream>
#include <mutex>
#include <sstream>

// Session owns shared annotation/configuration; workers own per-read scratch.

namespace nustar::star2711b::solo {
using detail::AnnotationIndex;
using detail::Context;
using detail::Counter;
namespace FeatureTypes = detail::FeatureTypes;
using detail::UnsignedWord;
using detail::WorkerState;
struct Session::Impl {
    Context parameters;
    std::unique_ptr<AnnotationIndex> transcriptome;
    std::vector<std::shared_ptr<Worker::Impl>> workers;
    std::size_t next_worker = 0;
    std::vector<detail::ReadCorrection> corrections;
    std::mutex mutex;
    std::string prefix;
    const ReferenceMetadata& metadata;
    Impl(const Config& config, const std::filesystem::path& genome, std::string out,
         unsigned threads, unsigned mates, const ReferenceMetadata& ref)
        : prefix(std::move(out)), metadata(ref) {
        auto& p = parameters;
        p.inOut->logStdOut = &std::cout;
        p.inOut->logMain.open(prefix + "Log.out");
        p.readOptions(config.arguments);
        p.runMode = "alignReads";
        p.runThreadN = threads;
        p.runDirPerm = 0700;
        p.runRestart.type = 0;
        p.quant = {};
        p.outSAMattrPresent = {};
        p.outBAMcoord = config.coordinate_bam;
        p.readFilesTypeN = config.bam_input ? 10 : 2;
        p.outputAttributes = config.output_attributes;
        for (const auto& tag : config.output_attributes) {
            if (tag == "CR")
                p.outSAMattrPresent.CR = true;
            if (tag == "CY")
                p.outSAMattrPresent.CY = true;
            if (tag == "UR")
                p.outSAMattrPresent.UR = true;
            if (tag == "UY")
                p.outSAMattrPresent.UY = true;
            if (tag == "CB")
                p.outSAMattrPresent.CB = true;
            if (tag == "UB")
                p.outSAMattrPresent.UB = true;
            if (tag == "sM")
                p.outSAMattrPresent.sM = true;
            if (tag == "sS")
                p.outSAMattrPresent.sS = true;
            if (tag == "sQ")
                p.outSAMattrPresent.sQ = true;
            if (tag == "sF")
                p.outSAMattrPresent.sF = true;
        }
        p.readNmates = mates;
        p.readNends = mates + (config.type == "SmartSeq" || config.barcode_mate ? 0 : 1);
        for (unsigned end = 0; end < 2; ++end) {
            p.pClip.in[end].N.resize(mates);
            for (unsigned m = 0; m < mates; ++m)
                p.pClip.in[end].N[m] = config.clipping[m][end];
        }
        p.pGe.gDir = genome.string();
        p.pGe.transform.outQuant = false;
        p.outFileNamePrefix = prefix;
        p.outFileTmp = prefix + "_STARtmp/";
        std::filesystem::remove_all(p.outFileTmp);
        std::filesystem::create_directories(p.outFileTmp);
        for (std::size_t i = 0; i < ref.chromosomes.size(); ++i)
            if (ref.chromosomes[i].name == "chrM" || ref.chromosomes[i].name == "MT")
                p.pGe.chrSet.mito.insert(i);

        // SmartSeq CBs are read-file identifiers. One input group uses one CB.
        p.outSAMattrRG = config.read_groups;
        p.options.initialize(&p);
        transcriptome = std::make_unique<AnnotationIndex>(p);
        // STAR constructs every SoloRead before mapping begins. In complex
        // barcode mode the first read determines umiL, so lazy worker creation
        // would initialize later workers with different homopolymer constants.
        for (unsigned id = 0; id < threads; ++id)
            workers.push_back(std::make_shared<Worker::Impl>(p, *transcriptome, id));
    }
    ~Impl() {
        workers.clear();
        transcriptome.reset();
        parameters.ownedStreams.clear();
        std::error_code ignored;
        std::filesystem::remove_all(parameters.outFileTmp, ignored);
    }
};
Session::Session(const Config& c, const std::filesystem::path& g, std::string prefix,
                 unsigned threads, unsigned mates, const ReferenceMetadata& ref)
    : impl_(std::make_unique<Impl>(c, g, std::move(prefix), threads, mates, ref)) {}
Session::~Session() = default;
std::unique_ptr<Worker> Session::createWorker() {
    std::lock_guard lock(impl_->mutex);
    return std::make_unique<Worker>(impl_->workers[impl_->next_worker++]);
}
void Session::finish() {
    auto& state = *impl_;
    auto& p = state.parameters;
    std::vector<WorkerState*> workers;
    p.countingStats = {};
    for (auto& worker : state.workers) {
        workers.push_back(&worker->native);
        p.countingStats.readN += worker->statistics.reads;
        p.countingStats.mappedReadsU += worker->statistics.mapped_unique;
        p.countingStats.mappedReadsM += worker->statistics.mapped_multiple;
    }
    // outputSJ uses filtered junctions for the Solo SJ feature.
    // The upstream counter also opens this file when the junction list is empty.
    std::ofstream sj_temporary;
    if (p.options.featureYes[FeatureTypes::SJ])
        sj_temporary.open(p.outFileTmp + "SJ.start_gap.tsv");
    std::ifstream sj(state.prefix + "SJ.out.tab");
    std::string line, chr;
    while (std::getline(sj, line)) {
        std::istringstream row(line);
        UnsignedWord start, end;
        if (!(row >> chr >> start >> end))
            continue;
        for (const auto& c : state.metadata.chromosomes)
            if (c.name == chr) {
                p.sjAll[0].push_back(c.start.value + start - 1);
                p.sjAll[1].push_back(end - start + 1);
                if (sj_temporary.is_open())
                    sj_temporary << p.sjAll[0].back() << '\t' << p.sjAll[1].back() << '\n';
                break;
            }
    }
    if (sj_temporary.is_open())
        sj_temporary.close();
    Counter counter(workers.data(), p, *state.transcriptome);
    counter.processAndOutput();
    if (p.options.samAttrYes)
        state.corrections =
            std::move(counter.features_[p.options.featureInd[p.options.samAttrFeature]]->readInfo);
}
bool Session::hasCorrections() const { return impl_->parameters.options.samAttrYes; }
std::pair<std::string, std::string> Session::correctedTags(std::uint64_t read) const {
    const auto& state = *impl_;
    if (read >= state.corrections.size())
        return {"-", "-"};
    const auto& c = state.corrections[read];
    const auto& options = state.parameters.options;
    return {c.cb == detail::UnsignedWord(-1) ? "-" : options.cbWLstr[c.cb],
            c.umi == std::uint32_t(-1) ? "-"
                                       : detail::convertNuclInt64toString(c.umi, options.umiL)};
}
} // namespace nustar::star2711b::solo
