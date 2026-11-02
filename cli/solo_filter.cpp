#include "solo_filter.hpp"
#include "../solo/session.hpp"
#include <iostream>
#include <stdexcept>

int runSoloCellFiltering(int argc,char** argv) {
    using namespace nustar::star2711b;
    try {
        solo::Config config;
        std::filesystem::path input;
        std::string prefix,log_prefix="./";
        for(int i=1;i<argc;++i) {
            const std::string flag=argv[i];
            if(flag=="--runMode" && i+3<argc) {
                if(std::string_view(argv[++i])!="soloCellFiltering") throw std::invalid_argument("invalid runMode");
                input=argv[++i];prefix=argv[++i];
            } else if(flag=="--outFileNamePrefix" && i+1<argc) log_prefix=argv[++i];
            else if(flag.starts_with("--solo") || flag=="--runThreadN" || flag=="--runRNGseed") {
                config.arguments.push_back(flag);
                if(i+1==argc || std::string_view(argv[i+1]).starts_with("--")) throw std::invalid_argument("missing value: "+flag);
                while(i+1<argc && !std::string_view(argv[i+1]).starts_with("--")) config.arguments.emplace_back(argv[++i]);
            } else throw std::invalid_argument("unsupported soloCellFiltering option: "+flag);
        }
        if(input.empty() || prefix.empty()) throw std::invalid_argument("--runMode soloCellFiltering requires input directory and output prefix");
        auto result=solo::filterCells(config,input,prefix,log_prefix);
        if(auto failure=std::get_if<Failure>(&result)) throw std::runtime_error(failure->message);
        return 0;
    } catch(const std::exception& error) {
        std::cerr<<"soloCellFiltering: "<<error.what()<<'\n';return 2;
    }
}
