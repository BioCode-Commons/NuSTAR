#pragma once

#include <string>
#include <time.h>

namespace nustar::star2711b::solo::detail {

std::string timeMonthDayTime();
std::string timeMonthDayTime(time_t& rawTime);

} // namespace nustar::star2711b::solo::detail
