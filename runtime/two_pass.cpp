#include "two_pass.hpp"
#include "../index/junction_insert.hpp"
#include "../annotation/junction_database.hpp"
#include "../annotation/table_writer.hpp"
#include "../reference/synthetic_junction.hpp"
#include <algorithm>
#include <fstream>
#include <chrono>
#include <iostream>
#include <sstream>
#include <unordered_map>

namespace nustar::star2711b {
namespace {
template<class T>T checked(Result<T> r) {
    if(auto f=std::get_if<Failure>(&r)) throw std::runtime_error(f->message);
    return std::get<T>(std::move(r));
}
AlignmentResources load(const std::filesystem::path& dir) {
    auto ref=checked(loadGenomeReference(dir));
    SuffixArrayIndexBackend backend;
    auto index=checked(backend.load({"nustar-sa",dir/"index.sa",{}},ref));
    return {std::move(ref),std::move(index)};
}
AlignmentResources insert(const AlignmentResources& old, const RunControl& control,
    const std::filesystem::path& junctions, const std::filesystem::path& out, unsigned threads) {
    const auto& reference=*old.reference;
    const auto& chromosomes=reference.storage().metadata.chromosomes;
    std::unordered_map<std::string,std::size_t> by_name;
    for(std::size_t i=0;i<chromosomes.size();++i) by_name.emplace(chromosomes[i].name,i);
    std::vector<JunctionProposal> proposals;
    for(const auto& j:reference.junctions()) {
        auto chr=std::upper_bound(chromosomes.begin(),chromosomes.end(),j.reported_intron.begin.value,
            [](auto pos,const auto& c){return pos<c.start.value;})-chromosomes.begin()-1;
        proposals.push_back({ChromosomeId{static_cast<std::uint32_t>(chr)},j.reported_intron,j.strand,
            JunctionOrigin::ExistingIndex,30,proposals.size(),{}});
    }
    std::ifstream file(junctions);
    if(!file) throw std::runtime_error("cannot open pass-1 junctions: "+junctions.string());
    std::string line,chr;
    while(std::getline(file,line)) {
        std::istringstream row(line);
        std::uint64_t start,end; unsigned strand;
        if(!(row>>chr>>start>>end>>strand)) continue;
        const auto c=by_name.at(chr); const auto origin=chromosomes[c].start.value;
        proposals.push_back({ChromosomeId{static_cast<std::uint32_t>(c)},
            {GenomeOffset{origin+start-1},GenomeOffset{origin+end}},Strand(strand),
            JunctionOrigin::JunctionFile,0,proposals.size(),{}});
    }
    AnnotationConfig config{};
    config.sjdbOverhang=reference.junctions().empty() ? control.sjdb_overhang : reference.sjdbOverhang();
    config.limitSjdbInsertNsj=UINT64_MAX;
    GenomeJunctionPreparer preparer;
    auto prepared=checked(preparer.prepare(proposals,reference.view(),config));
    if(prepared.entries.size()>control.junction_limit)
        throw std::runtime_error("junction count exceeds limitSjdbInsertNsj");
    auto synthetic=makeSyntheticJunctionSequence(prepared.entries,reference.storage().forward,config.sjdbOverhang);
    auto updated=std::make_shared<const GenomeReferenceData>(reference.storageLease(),
        std::move(prepared.entries),std::move(synthetic),config.sjdbOverhang,!control.save_inserted_index);
    const auto& index=dynamic_cast<const SuffixArrayIndex&>(*old.index);
    auto new_index=checked(insertJunctions(index,reference,updated,threads));
    std::filesystem::create_directories(out);
    checked(writeJunctionDatabaseTables(updated->view(),out));
    if(control.save_inserted_index) {
        checked(saveGenomeReference(*updated,out));
        SuffixArrayIndexBackend backend;
        checked(backend.save(*new_index,out/"index.sa"));
        const auto base=index.artifact().parent_path();
        for(const auto* name:{"geneInfo.tab","transcriptInfo.tab","exonInfo.tab","exonGeTrInfo.tab","sjdbList.fromGTF.out.tab"})
            if(std::filesystem::exists(base/name))
                std::filesystem::copy_file(base/name,out/name,std::filesystem::copy_options::overwrite_existing);
    }
    return {std::move(updated),std::move(new_index)};
}
template<class Options,class Run>
auto run(const Options& options,Run align) -> decltype(align(options)) {
    try {
        auto resources=options.run.resources.reference ? options.run.resources : load(options.genome_dir);
        const auto prefix=options.outFileNamePrefix.value_or(options.output_sam.string()+".");
        const auto directory=prefix+"_STARpass1/";
        std::filesystem::create_directories(directory);
        auto first=options;
        first.run.resources=resources; first.run.twopass=false;
        first.run.sam_enabled=false;
        first.run.bam.unsorted=false;
        first.run.bam.coordinate=false;
        first.run.solo_config.reset();
        first.run.read_limit=std::min(options.run.pass1_reads,options.run.read_limit);
        first.outFileNamePrefix=directory;
        first.output_sam=directory+"Aligned.out.sam";
        first.outSAMunmappedWithin=false;
        const auto pass1_start=std::chrono::steady_clock::now();
        auto result=align(first);
        const auto insertion_start=std::chrono::steady_clock::now();
        if(auto failure=std::get_if<Failure>(&result)) return *failure;
        auto second=options;
        second.run.twopass=false;
        second.run.resources=insert(resources,options.run,directory+"SJ.out.tab",prefix+"_STARgenome/",options.threads);
        first.run.resources={}; resources={};
        const auto pass2_start=std::chrono::steady_clock::now();
        auto final=align(second);
        const auto finished=std::chrono::steady_clock::now();
        std::clog<<"twopass: pass1_seconds="<<std::chrono::duration<double>(insertion_start-pass1_start).count()
            <<" junction_insert_seconds="<<std::chrono::duration<double>(pass2_start-insertion_start).count()
            <<" pass2_seconds="<<std::chrono::duration<double>(finished-pass2_start).count()<<'\n';
        return final;
    } catch(const std::exception& e) {return Failure{ErrorKind::Index,Stage::Sjdb,{},1,e.what()};}
}
}
Result<SingleSamReport> runTwoPass(const SingleSamOptions& options) {return run(options,runSingleSam);}
Result<PairedSamReport> runTwoPass(const PairedSamOptions& options) {return run(options,runPairedSam);}
}
