#pragma once

#include "../configuration/context.hpp"
#include <fstream>

namespace nustar::star2711b::solo::detail {

std::fstream& fstrOpen(std::string fileName, std::string errorID, Context& context_,
                       bool flagDelete);
std::ofstream& ofstrOpen(std::string fileName, std::string errorID, Context& context_);
std::ifstream& ifstrOpen(std::string fileName, std::string errorID, std::string solutionString,
                         Context& context_);
std::ifstream& ifstrOpenGenomeFile(std::string fileName, std::string errorID, Context& context_);
void createDirectory(const std::string dirPathIn, const mode_t dirPerm,
                     const std::string dirParameter, Context& context_);

} // namespace nustar::star2711b::solo::detail
