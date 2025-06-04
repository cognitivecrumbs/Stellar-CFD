#include <ThreadPool.h>

ThreadPool::ThreadPool(size_t numThreads) : stop(false) {
    start(numThreads);
}

ThreadPool::~ThreadPool() {
    if (!stop) {
        stopPool();
    }
}

void ThreadPool::restart(size_t numThreads) {
    stopPool();
    start(numThreads);
}

void ThreadPool::stopPool() {
    {
        std::unique_lock<std::mutex> lock(queueMutex);
        stop = true;
    }
    condition.notify_all();
    for (std::thread& worker : workers) {
        worker.join();
    }
    workers.clear();
}

void ThreadPool::start(size_t numThreads) {
    stop = false;
    workers.reserve(numThreads);
    for (size_t i = 0; i < numThreads; ++i) {
        workers.emplace_back([this] {
            while (true) {
                std::function<void()> task;
                {
                    std::unique_lock<std::mutex> lock(this->queueMutex);
                    this->condition.wait(lock, [this] { return this->stop || !this->taskQueue.empty(); });
                    if (this->stop && this->taskQueue.empty()) {
                        return;
                    }
                    task = std::move(this->taskQueue.front());
                    this->taskQueue.pop();
                }
                task();
            }
        });
    }
}

// template <class F>
// void ThreadPool::enqueue(F&& f) {
//     {
//         std::unique_lock<std::mutex> lock(queueMutex);
//         if (stop) {
//             throw std::runtime_error("ThreadPool is stopped, cannot enqueue tasks");
//         }
//         taskQueue.push(std::forward<F>(f));
//     }
//     condition.notify_one();
// }