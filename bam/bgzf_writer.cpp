#include "bgzf_writer.hpp"
#include <zlib.h>

namespace nustar::star2711b::bam {
void BgzfWriter::compress(std::span<const std::byte> input, Buffer& output, int level) {
    struct Deflater {
        z_stream stream{};
        int level;
        explicit Deflater(int value) : level(value) {
            if (deflateInit2(&stream, level, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) != Z_OK)
                throw std::runtime_error("cannot initialize BAM compression");
        }
        ~Deflater() { deflateEnd(&stream); }
    };
        thread_local std::unique_ptr<Deflater> codec;
        if (!codec || codec->level != level) codec = std::make_unique<Deflater>(level);
        auto& stream = codec->stream;
        if (deflateReset(&stream) != Z_OK) throw std::runtime_error("cannot reset BAM compressor");
        const auto begin = output.size;
        if (begin + 65536 > output.capacity)
            output.reserve(std::max(begin + 65536, output.capacity <= SIZE_MAX/2 ? output.capacity*2 : begin+65536));
        auto* dst = reinterpret_cast<unsigned char*>(output.data.get() + begin);
        constexpr unsigned char header[18] = {31,139,8,4,0,0,0,0,0,255,6,0,66,67,2,0,0,0};
        std::memcpy(dst, header, sizeof(header));
        stream.next_in = reinterpret_cast<Bytef*>(const_cast<std::byte*>(input.data()));
        stream.avail_in = input.size(); stream.next_out = dst + 18; stream.avail_out = 65536 - 26;
        if (deflate(&stream, Z_FINISH) != Z_STREAM_END)
            throw std::runtime_error("compressed BAM block exceeds BGZF size limit");
        const std::size_t n = stream.total_out + 26;
        const auto put = [](unsigned char* p, std::uint32_t value, unsigned bytes) {
            for (unsigned i = 0; i < bytes; ++i) p[i] = value >> (8*i);
        };
        put(dst + 16, n - 1, 2);
        put(dst + n - 8, crc32(0, reinterpret_cast<const Bytef*>(input.data()), input.size()), 4);
        put(dst + n - 4, input.size(), 4);
        output.size = begin + n;
    }
}
