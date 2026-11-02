#include "error.hpp"
#include "time.hpp"
#include <mutex>

namespace nustar::star2711b::solo::detail {

/*
    functions that handle errors and warnings
*/

void exitWithError(std::string messageOut, std::ostream& streamOut1, std::ostream& streamOut2,
                   int errorInt, Context& context_) {
    std::lock_guard lock(context_.errorMutex);
    time_t timeCurrent;
    time(&timeCurrent);
    if (streamOut1.good()) {
        streamOut1 << "\n"
                   << messageOut << std::endl
                   << timeMonthDayTime(timeCurrent) << " ...... FATAL ERROR, exiting\n"
                   << std::flush;
    };
    if (streamOut2.good()) {
        streamOut2 << "\n"
                   << messageOut << std::endl
                   << timeMonthDayTime(timeCurrent) << " ...... FATAL ERROR, exiting\n"
                   << std::flush;
    };
    context_.inOut.reset(); // close logs before process termination
    exit(errorInt);
};

void warningMessage(std::string messageOut, std::ostream& streamOut1, std::ostream& streamOut2,
                    Context& context_) {
    std::lock_guard lock(context_.errorMutex);
    if (streamOut1.good()) {
        streamOut1 << "!!!!! WARNING: " << messageOut << std::endl;
    };
    if (streamOut2.good()) {
        streamOut2 << "!!!!! WARNING: " << messageOut << std::endl;
    };
};

} // namespace nustar::star2711b::solo::detail
