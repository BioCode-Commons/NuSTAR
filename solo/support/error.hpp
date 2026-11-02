#pragma once

#include "../configuration/context.hpp"
#include "types.hpp"

namespace nustar::star2711b::solo::detail {

void exitWithError(std::string messageOut, std::ostream& streamOut1, std::ostream& streamOut2,
                   int errorInt, Context& context_);
void warningMessage(std::string messageOut, std::ostream& streamOut1, std::ostream& streamOut2,
                    Context& context_);

} // namespace nustar::star2711b::solo::detail
