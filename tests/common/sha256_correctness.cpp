#include "../../common/sha256.hpp"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace nustar::star2711b;
int main() {
    try {
        // Independent hashlib SHA-256 known answers, including both padding
        // branches, unaligned input and many complete compression blocks.
        const std::pair<std::size_t, const char*> cases[] = {
            {0, "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"},
            {1, "e7cf46a078fed4fafd0b5e3aff144802b853f8ae459a4f0c14add3314b7cc3a6"},
            {55, "2900465fcb533e05a158fd2b3be0e5e3b03740d83060aa3580e0d98a96bf2384"},
            {56, "31454ff48ef36af2f08fd511bdc37d9d5855ac23e992e5ff5445cb6b7674a674"},
            {63, "5f6401b96532c36de4e65beec0409b69b1d181864c8009b7a04f43e5d56350d1"},
            {64, "94eb5de4943613fd048dc93393ab06877405faa39c11f53e9386083339833e7e"},
            {65, "fc518669b6eb4b4dd91827ecacef86689c725bd5bab888fd3b26dbb196eec954"},
            {119, "b0dc41b1a384e2f1203f0351b38fbeaafceef577ce1191d5bfc25da39f721eae"},
            {120, "5df24dd802ac26132ce608dcb5f09841eef039ee0f152acf98d26d17fe4e88e6"},
            {127, "0fe729ff19257bd6fec853acc2ea355f6b34b58e6c0f684c3e188fcdfcd9baae"},
            {128, "0aedd4856f8eba0963627336ad5144a9a7dbe12498e6066f0165fc97d8ddee4c"},
            {129, "4f1757ae4bffbae86d775b831765b75af154d52f7deaa46dd378051a2d3ad57f"},
            {1024, "ffbad8f947474cfdd5b2bb22d7e0bf5ee8ba2b7af859d0c2bb28622db6a4be47"},
            {1000000, "047207f05c235c6490b408242a86a6f1d4c03d0e754902bd13efeb29f7f08a66"},
        };
        for (const auto& [length, expected] : cases) {
            std::vector<std::byte> data(length+1);
            for (std::size_t i = 0; i < length; ++i) data[i+1] = std::byte((i*37+11)%256);
            const auto bytes = std::span(data).subspan(1);
            if (sha256(bytes) != expected) throw std::runtime_error("whole-buffer digest mismatch");
            for (const std::size_t chunk : {1, 7, 63, 64, 65, 1024}) {
                Sha256 hash;
                for (std::size_t offset = 0; offset < length; offset += chunk)
                    hash.update(bytes.subspan(offset, std::min(chunk, length-offset)));
                if (hash.finish() != expected) throw std::runtime_error("streaming digest mismatch");
            }
        }
        return 0;
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
