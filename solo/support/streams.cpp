#include "streams.hpp"
#include "error.hpp"
#include <cstring>
#include <fstream>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/statvfs.h>

namespace nustar::star2711b::solo::detail {

void createDirectory(const std::string dirPathIn, const mode_t dirPerm,
                     const std::string dirParameter, Context& context_) {
    std::string dirPath = dirPathIn.substr(0, dirPathIn.find_last_of('/') + 1);
    if (mkdir(dirPath.c_str(), dirPerm) == -1) {
        if (errno == EEXIST) { // directory exists
            context_.inOut->logMain << dirParameter
                                    << " directory exists and will be overwritten: " << dirPath
                                    << std::endl;
        } else { // other error
            // will try to create parent directories
            size_t i1 = dirPath.find_first_of('/', 1);
            while (i1 < dirPath.size()) {
                std::string dirPath1 = dirPath.substr(0, i1);
                if (mkdir(dirPath1.c_str(), dirPerm) == -1) {
                    if (!(errno == EEXIST)) { // error
                        exitWithError("EXITING because of fatal OUTPUT FILE error: could not "
                                      "create output directory: " +
                                          dirPath1 + " for " + dirParameter + " " + dirPathIn +
                                          "\n ERROR: " + strerror(errno) +
                                          "\nSOLUTION: check the path and permissions.\n",
                                      std::cerr, context_.inOut->logMain, EXIT_CODE_PARAMETER,
                                      context_);
                    };
                };
                i1 = dirPath.find_first_of('/', i1 + 1);
            };
            context_.inOut->logMain
                << dirParameter << " directory and its parents created: " << dirPath << std::endl;
        };
    } else {
        context_.inOut->logMain << dirParameter << " directory created: " << dirPath << std::endl;
    };
};

std::ofstream& ofstrOpen(std::string fileName, std::string errorID,
                         Context& context_) { // open file 'fileName', generate error if cannot open
    std::ofstream& ofStream =
        *new std::ofstream(fileName.c_str(), std::fstream::out | std::fstream::trunc);
    if (ofStream.fail()) { //
        std::ostringstream errOut;
        errOut << errorID
               << ": exiting because of *OUTPUT FILE* error: could not create output file "
               << fileName << "\n";
        errOut << "SOLUTION: check that the path exists and you have write permission for this "
                  "file. Also check "
                  "ulimit -n"
                  " and increase it to allow more open files.\n";
        exitWithError(errOut.str(), std::cerr, context_.inOut->logMain, EXIT_CODE_FILE_OPEN,
                      context_);
    };
    {
        std::lock_guard<std::mutex> lock(context_.streamMutex);
        context_.ownedStreams.emplace_back(&ofStream);
    }
    return ofStream;
};

std::fstream& fstrOpen(std::string fileName, std::string errorID, Context& context_,
                       bool flagDelete) { // open file 'fileName', generate error if cannot open

    std::fstream* fStreamP;
    if (flagDelete) { // truncate the file if it exists
        fStreamP = new std::fstream(fileName.c_str(),
                                    std::fstream::in | std::fstream::out | std::fstream::trunc);
    } else { // try to open exising file
        fStreamP = new std::fstream(fileName.c_str(), std::fstream::in | std::fstream::out);
        if (fStreamP->fail()) {
            fStreamP->clear();
            fStreamP->open(fileName.c_str(),
                           std::fstream::in | std::fstream::out | std::fstream::trunc);
        }
    };

    if (fStreamP->fail()) { //
        std::ostringstream errOut;
        errOut << errorID
               << ": exiting because of *OUTPUT FILE* error: could not create input/output file "
               << fileName << "\n";
        errOut
            << "Solution: check that the path exists and you have write permission for this file\n";
        exitWithError(errOut.str(), std::cerr, context_.inOut->logMain, EXIT_CODE_FILE_OPEN,
                      context_);
    };
    {
        std::lock_guard<std::mutex> lock(context_.streamMutex);
        context_.ownedStreams.emplace_back(fStreamP);
    }
    return *fStreamP;
};

std::ifstream& ifstrOpen(std::string fileName, std::string errorID, std::string solutionString,
                         Context& context_) {
    // open file 'fileName', generate error if cannot open
    std::ifstream& ifStream = *new std::ifstream(fileName.c_str());
    if (ifStream.fail()) { //
        std::ostringstream errOut;
        errOut << errorID << ": exiting because of *INPUT FILE* error: could not open input file "
               << fileName << "\n";
        errOut
            << "Solution: check that the file exists and you have read permission for this file\n";
        if (solutionString.size() > 0) {
            errOut << "          " << solutionString << "\n";
        };
        exitWithError(errOut.str(), std::cerr, context_.inOut->logMain, EXIT_CODE_FILE_OPEN,
                      context_);
    };
    {
        std::lock_guard<std::mutex> lock(context_.streamMutex);
        context_.ownedStreams.emplace_back(&ifStream);
    }
    return ifStream;
};

std::ifstream& ifstrOpenGenomeFile(std::string fileName, std::string errorID, Context& context_) {
    // open one of the genome files
    return ifstrOpen(context_.pGe.gDir + "/" + fileName, errorID,
                     "if this file is missing from the genome directory, you will need to "
                     "*re-generate the genome*",
                     context_);
};

void copyFile(std::string fileIn, std::string fileOut) { // copy fileIn into FileOut
    std::ifstream src(fileIn, std::ios::binary);
    std::ofstream dst(fileOut, std::ios::binary);
    dst << src.rdbuf();
};

} // namespace nustar::star2711b::solo::detail
