#include "../../bam/sort.hpp"
#include <algorithm>
#include <atomic>
#include <cstring>
#include <iostream>
#include <random>
#include <thread>
#include <zlib.h>

using namespace nustar::star2711b;
using namespace nustar::star2711b::bam;
void require(bool yes, const char* message) { if (!yes) throw std::runtime_error(message); }
Bytes record(std::uint32_t tid, std::uint32_t pos, std::uint64_t marker, unsigned length = 96) {
    Bytes bytes;
    append32(bytes, length - 4); append32(bytes, tid); append32(bytes, pos);
    bytes.resize(length);
    std::memcpy(bytes.data() + 36, &marker, 8);
    return bytes;
}
struct Expected {
    std::uint32_t tid, pos;
    std::uint64_t read, order, marker;
    auto key() const { return std::tuple{tid, pos, read, order}; }
};
void scenario(const std::filesystem::path& root, unsigned threads, unsigned mode, std::uint64_t budget) {
    const auto directory = root / (std::to_string(threads) + "-" + std::to_string(mode));
    std::filesystem::create_directories(directory);
    const unsigned count = mode == 0 ? 0 : mode == 1 ? 1 : mode == 2 ? 20000 : 180000;
    std::vector<OutputBatch> batches;
    std::vector<Expected> expected;
    std::mt19937 random(731);
    for (unsigned i = 0; i < count; ++i) {
        if (i % 1024 == 0) batches.push_back({BatchId{batches.size()}, {}});
        const std::uint64_t read = mode == 4 ? 0 :
            ((mode == 2 || mode == 3) && (i & 1) ? 1ULL << 63 : 0) + count - i;
        EncodedFragment fragment{{}, {}, {}}; fragment.fragment.ordinal = read;
        for (unsigned order = 0; order < (mode == 4 ? 1U : 2U); ++order) {
            const std::uint32_t tid = mode == 4 ? UINT32_MAX : mode == 3 && i > count / 4 ? 100U : random() % 5;
            const std::uint32_t pos = mode == 4 ? UINT32_MAX : mode == 3 ? 7U : random() % 11;
            const auto marker = mode == 4 ? std::uint64_t(i) : 2ULL * i + order;
            auto bytes = record(tid, pos, marker, mode == 5 && i == 99 ? 2000000 : 96);
            fragment.sam.insert(fragment.sam.end(), bytes.begin(), bytes.end());
            expected.push_back({tid, pos, read, order, marker});
        }
        batches.back().fragments.push_back(std::move(fragment));
    }
    std::stable_sort(expected.begin(), expected.end(), [](const auto& a, const auto& b) { return a.key() < b.key(); });
    {
        CoordinateSorter sorter(directory, budget, threads);
        std::atomic<std::size_t> next{0};
        std::exception_ptr failure;
        std::mutex mutex;
        std::vector<std::jthread> workers;
        for (unsigned i = 0; i < threads; ++i) workers.emplace_back([&, i] {
            // Distinct payloads with completely equal keys must retain their
            // submitted order, including queued versus direct spool writes.
            if (mode == 4 && i) return;
            try {
                for (;;) {
                    auto j = next.fetch_add(1);
                    if (j >= batches.size()) break;
                    sorter.submit(batches[j]);
                }
            } catch (...) { std::lock_guard lock(mutex); failure = std::current_exception(); }
        });
        workers.clear();
        if (failure) std::rethrow_exception(failure);
        sorter.finish();
        std::size_t at = 0;
        const auto verify = [&](std::uint64_t read, std::span<const std::byte> bytes) {
            require(at < expected.size(), "extra record");
            const auto& e = expected[at++];
            std::uint64_t marker; std::memcpy(&marker, bytes.data() + 36, 8);
            require(read == e.read && read32(bytes.data() + 4) == e.tid &&
                    read32(bytes.data() + 8) == e.pos && marker == e.marker, "sort differs from reference");
        };
        if (mode == 5) {
            const auto path = directory / "long-records.bgzf";
            sorter.writeBam(path.string(), {}, 1);
            gzFile file = gzopen(path.c_str(), "rb");
            require(file != nullptr, "cannot open sorted BGZF");
            std::array<std::byte, 4> length;
            while (gzread(file, length.data(), length.size()) == 4) {
                const auto n = std::uint64_t(read32(length.data())) + 4;
                Bytes record(n); std::copy(length.begin(), length.end(), record.begin());
                require(gzread(file, record.data()+4, n-4) == std::int64_t(n-4), "truncated sorted BGZF");
                require(at < expected.size(), "extra sorted BGZF record");
                verify(expected[at].read, record);
            }
            require(gzeof(file) && gzclose(file) == Z_OK, "BGZF decoding failed");
        } else if (mode == 2 || mode == 4) {
            Bytes pending;
            sorter.emit({}, [&](std::span<const std::byte> bytes) {
                pending.insert(pending.end(), bytes.begin(), bytes.end());
                std::size_t offset = 0;
                while (pending.size() - offset >= 4) {
                    const auto n = std::uint64_t(read32(pending.data() + offset)) + 4;
                    if (n > pending.size() - offset) break;
                    require(at < expected.size(), "extra bulk record");
                    verify(expected[at].read, std::span(pending).subspan(offset, n));
                    offset += n;
                }
                pending.erase(pending.begin(), pending.begin() + offset);
            });
            require(pending.empty(), "truncated bulk record");
        } else sorter.emit(verify);
        require(at == expected.size(), "lost record");
        auto metrics = sorter.metrics();
        require(metrics.records == expected.size(), "incorrect record count");
        if (mode >= 3) require(metrics.spilled_bytes > 0, "spill path not exercised");
        if (mode == 2) require(metrics.spilled_bytes == 0, "unexpected in-memory spill");
        std::cout << "PASS mode=" << mode << " threads=" << threads << " records=" << at
                  << " written=" << metrics.spilled_bytes << " partitions=" << metrics.partitions
                  << " radix_s=" << metrics.radix_seconds << " pack_s=" << metrics.pack_seconds << '\n';
        std::cout << "IO read_calls=" << metrics.temporary_reads.calls
                  << " write_calls=" << metrics.temporary_writes.calls
                  << " read_bytes=" << metrics.read_bytes
                  << " output_calls=" << metrics.output_writes.calls
                  << " read_worker_s=" << metrics.read_seconds
                  << " parse_worker_s=" << metrics.parse_seconds
                  << " compression_worker_s=" << metrics.compression_seconds
                  << " wait_wall_s=" << metrics.wait_seconds
                  << " sync_wall_s=" << metrics.sync_seconds << '\n';
    }
    std::filesystem::remove_all(directory);
}
int main(int argc, char** argv) {
    try {
        require(argc == 2, "expected scratch directory");
        const std::filesystem::path root = argv[1];
        for (unsigned threads : {1U, 2U, 4U, 32U})
            for (unsigned mode = 0; mode < 6; ++mode)
                scenario(root, threads, mode, mode < 3 ? 64ULL << 20 : 8ULL << 20);
        std::filesystem::create_directories(root / "errors");
        {
            CoordinateSorter sorter(root / "errors", 8ULL << 20, 4);
            OutputBatch batch{{}, {}};
            batch.fragments.push_back({{}, Bytes(7), {}});
            bool failed = false;
            try { sorter.submit(batch); } catch (const std::exception&) { failed = true; }
            require(failed, "malformed record accepted");
            failed = false;
            try { sorter.finish(); } catch (const std::exception&) { failed = true; }
            require(failed, "submit failure not retained");
        }
        {
            CoordinateSorter sorter(root / "errors", 8ULL << 20, 4);
            sorter.cancel(); bool failed = false;
            try { sorter.finish(); } catch (const std::exception&) { failed = true; }
            require(failed, "cancel ignored");
        }
        {
            CoordinateSorter sorter(root / "errors", 8ULL << 20, 4);
            OutputBatch batch{{}, {{{}, record(0, 0, 0), {}}}};
            sorter.submit(batch);
            bool failed = false;
            try { sorter.emit([](auto, auto) { throw std::runtime_error("injected output failure"); }); }
            catch (const std::exception&) { failed = true; }
            require(failed, "output failure ignored");
        }
        {
            CoordinateSorter sorter(root / "missing" / "directory", 8ULL << 20, 4);
            OutputBatch batch{{}, {{{}, record(0, 0, 0, 6000000), {}}}};
            bool failed = false;
            try { sorter.submit(batch); } catch (const std::exception&) { failed = true; }
            require(failed, "temporary write failure ignored");
        }
        {
            CoordinateSorter sorter(root / "errors", 8ULL << 20, 4);
            OutputBatch batch{{}, {{{}, record(0, 0, 0, 6000000), {}}}};
            sorter.submit(batch); sorter.finish();
            for (const auto& file : std::filesystem::directory_iterator(root / "errors"))
                std::filesystem::resize_file(file.path(), 1);
            bool failed = false;
            try { sorter.emit([](auto, auto) {}); } catch (const std::exception&) { failed = true; }
            require(failed, "truncated bucket ignored");
        }
        std::filesystem::remove_all(root / "errors");
        return 0;
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
