#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <source_location>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

namespace nustar::star2711b::solo::detail {

// STAR packs and compares records as unsigned long long. Preserve the exact
// type (including literal promotion), not merely its width.
using UnsignedWord = unsigned long long;
using SignedWord = long long;
static_assert(sizeof(UnsignedWord) == 8);
#ifdef COMPILE_FOR_LONG_READS
inline constexpr unsigned MAX_N_EXONS = 1000;
#else
inline constexpr unsigned MAX_N_EXONS = 20;
#endif
inline constexpr unsigned EXIT_CODE_PARAMETER = 102;
inline constexpr unsigned EXIT_CODE_INPUT_FILES = 104;
inline constexpr unsigned EXIT_CODE_GENOME_FILES = 105;
inline constexpr unsigned EXIT_CODE_FILE_OPEN = 109;
inline constexpr unsigned EXIT_CODE_FILE_WRITE = 110;
inline constexpr unsigned MARK_FRAG_SPACER_BASE = 11;
inline constexpr unsigned EX_R = 0;
inline constexpr unsigned EX_G = 1;
inline constexpr unsigned EX_L = 2;
inline constexpr unsigned EX_iFrag = 3;
inline constexpr unsigned EX_SIZE = 5;
inline std::string
sourceLocation(const std::source_location where = std::source_location::current()) {
    return std::string(where.file_name()) + ":" + std::to_string(where.line()) + ":" +
           where.function_name();
}

} // namespace nustar::star2711b::solo::detail
