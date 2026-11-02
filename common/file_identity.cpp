#include "file_identity.hpp"
#include "sha256.hpp"

#include <array>
#include <fstream>

namespace nustar::star2711b {

Result<FileIdentity> identifyFile(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) return Failure{ErrorKind::Io, Stage::Finish, std::nullopt, 1,
                               "cannot read file: " + path.string()};
    Sha256 hash;
    std::uint64_t bytes = 0;
    std::array<char, 1 << 16> buffer{};
    while (input) {
        input.read(buffer.data(), buffer.size());
        const auto count = input.gcount();
        if (count > 0) {
            hash.update(std::as_bytes(std::span(buffer.data(), static_cast<std::size_t>(count))));
            bytes += static_cast<std::uint64_t>(count);
        }
    }
    if (!input.eof()) return Failure{ErrorKind::Io, Stage::Finish, std::nullopt, 1,
                                    "error reading file: " + path.string()};
    return FileIdentity{path, bytes, hash.finish()};
}

} // namespace nustar::star2711b
