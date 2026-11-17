#pragma once
// Private BAM task pool. Only the coordinating thread waits for futures.
#include <condition_variable>
#include <deque>
#include <functional>
#include <future>
#include <mutex>
#include <thread>
#include <vector>

namespace nustar::star2711b::bam {
class WorkQueue {
    std::mutex mutex_;
    std::condition_variable ready_;
    std::deque<std::function<void()>> tasks_, urgent_;
    std::vector<std::thread> workers_;
    bool closed_ = false;
  public:
    explicit WorkQueue(unsigned threads) {
        try {
            for (unsigned i = 1; i < threads; ++i) workers_.emplace_back([this] {
                for (;;) {
                    std::function<void()> task;
                    {
                        std::unique_lock lock(mutex_);
                        ready_.wait(lock, [&] { return closed_ || !urgent_.empty() || !tasks_.empty(); });
                        auto& queue = urgent_.empty() ? tasks_ : urgent_;
                        if (queue.empty()) return;
                        task = std::move(queue.front()); queue.pop_front();
                    }
                    task();
                }
            });
        } catch (...) { close(); throw; }
    }
    ~WorkQueue() { close(); }
    void close() noexcept {
        { std::lock_guard lock(mutex_); closed_ = true; }
        ready_.notify_all();
        for (auto& worker : workers_) if (worker.joinable()) worker.join();
    }
    unsigned concurrency() const { return workers_.size() + 1; }
    template<class F> auto submit(F task, bool urgent = false) {
        using Result = std::invoke_result_t<F>;
        auto packaged = std::make_shared<std::packaged_task<Result()>>(std::move(task));
        auto future = packaged->get_future();
        if (workers_.empty()) (*packaged)();
        else {
            { std::lock_guard lock(mutex_);
              if (closed_) throw std::runtime_error("BAM task queue is closed");
              (urgent ? urgent_ : tasks_).push_back([packaged] { (*packaged)(); }); }
            ready_.notify_one();
        }
        return future;
    }
    template<class F> void parallel(unsigned count, F task) {
        std::vector<std::future<void>> futures;
        futures.reserve(count ? count-1 : 0);
        std::exception_ptr failure;
        try {
            for (unsigned i = 1; i < count; ++i) futures.push_back(submit([&, i] { task(i); }));
            task(0);
        } catch (...) { failure = std::current_exception(); }
        for (auto& future : futures) {
            try { future.get(); } catch (...) { if (!failure) failure = std::current_exception(); }
        }
        if (failure) std::rethrow_exception(failure);
    }
};
}
