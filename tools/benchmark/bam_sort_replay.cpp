// Replay identical deterministic BAM batches through old/new public sinks.
// Generation is included in wall time; submit_seconds is aggregate worker time.
#include "../../bam/io.hpp"
#include "../../bam/codec.hpp"
#include "../../third_party/htslib/htslib/hts.h"
#include <atomic>
#include <chrono>
#include <cstring>
#include <iostream>
#include <mutex>
#include <thread>

using namespace nustar::star2711b;
using Clock = std::chrono::steady_clock;
void checked(Status status) {
    if (auto* error = std::get_if<Failure>(&status)) throw std::runtime_error(error->message);
}
double seconds(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now()-start).count();
}
bam::Bytes record(std::uint64_t id) {
    bam::Bytes bytes;
    const std::uint32_t pos = (id * 2654435761ULL) % 1000000;
    const auto bin = std::uint32_t(hts_reg2bin(pos, pos + 150, 14, 5));
    for (std::uint32_t value : {269U, 0U, pos, (bin << 16) | (60U << 8) | 8U,
                               1U, 150U, UINT32_MAX, UINT32_MAX, 0U}) bam::append32(bytes, value);
    char name[8] = "r000000";
    auto number = id;
    for (int i = 6; i > 0; --i) { name[i] = '0'+number%10; number /= 10; }
    bytes.insert(bytes.end(), reinterpret_cast<std::byte*>(name), reinterpret_cast<std::byte*>(name)+8);
    bam::append32(bytes, 150U << 4);
    auto state = id + 123;
    const auto random = [&] { state ^= state << 13; state ^= state >> 7; state ^= state << 17; return state; };
    for (unsigned i = 0; i < 75; ++i) bytes.push_back(std::byte((1U << (random()%4)) | (1U << (4+random()%4))));
    for (unsigned i = 0; i < 150; ++i) bytes.push_back(std::byte(30 + random()%11));
    return bytes;
}
int main(int argc, char** argv) {
    try {
        if (argc != 5) throw std::runtime_error("usage: bam-sort-replay prefix records threads budget-bytes");
        const std::string prefix = argv[1];
        const auto count = std::stoull(argv[2]);
        const unsigned threads = std::stoul(argv[3]);
        if (!threads) throw std::runtime_error("threads must be positive");
        ReferenceMetadata reference{{{{0}, "chr1", {0}, 10000000}}, {10000000}, "test"};
        BamConfig config; config.coordinate = true; config.threads = threads;
        config.sort_memory = std::stoull(argv[4]);
        OutputFiles files; files.outFileNamePrefix = prefix; files.sam_enabled = false;
        files.job_started = files.mapping_started = std::chrono::system_clock::now();
        auto opened = openBamSink(files, config, reference);
        if (auto* error = std::get_if<Failure>(&opened)) throw std::runtime_error(error->message);
        auto sink = std::get<std::unique_ptr<IOutputSink>>(std::move(opened));
        const std::string header = "@HD\tVN:1.4\tSO:unsorted\n@SQ\tSN:chr1\tLN:10000000\n";
        checked(sink->begin({reinterpret_cast<const std::byte*>(header.data()), header.size()}));
        const auto start = Clock::now();
        std::atomic<std::uint64_t> next{0};
        std::mutex output_mutex, error_mutex;
        std::exception_ptr failure;
        std::vector<double> submit_seconds(threads);
        std::vector<std::jthread> workers;
        for (unsigned w = 0; w < threads; ++w) workers.emplace_back([&, w] {
            try {
                for (;;) {
                    const auto begin = next.fetch_add(4096);
                    if (begin >= count) break;
                    OutputBatch batch{BatchId{begin/4096}, {}};
                    for (auto i = begin; i < std::min<std::uint64_t>(count, begin+4096); ++i)
                        batch.fragments.push_back({FragmentId{0, i}, record(i), {}});
                    const auto submit = Clock::now();
                    if (sink->concurrentWrites()) checked(sink->submit(std::move(batch)));
                    else { std::lock_guard lock(output_mutex); checked(sink->submit(std::move(batch))); }
                    submit_seconds[w] += seconds(submit);
                }
            } catch (...) { std::lock_guard lock(error_mutex); failure = std::current_exception(); }
        });
        workers.clear();
        if (failure) std::rethrow_exception(failure);
        checked(sink->finish({}, {}));
        const auto collect = seconds(start);
        const auto finalize_start = Clock::now();
        checked(sink->finalize());
        double submit = 0;
        for (auto t : submit_seconds) submit += t;
        std::cout << "{\"records\":" << count << ",\"threads\":" << threads
                  << ",\"collect_wall_seconds\":" << collect << ",\"finalize_wall_seconds\":" << seconds(finalize_start)
                  << ",\"wall_seconds\":" << seconds(start) << ",\"submit_worker_seconds\":" << submit << "}\n";
        return 0;
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
