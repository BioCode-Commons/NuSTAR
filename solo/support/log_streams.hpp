#pragma once

#include <fstream>
#include <iostream>

namespace nustar::star2711b::solo::detail {

struct LogStreams {
    std::ostream* logStdOut = &std::cout;
    std::ofstream logMain;
};

} // namespace nustar::star2711b::solo::detail
