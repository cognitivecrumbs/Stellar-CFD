#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <vector>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <stdexcept>

class ThreadPool {
public:
    // Constructor and destructor
    ThreadPool(size_t numThreads);
    ~ThreadPool();

    // Template method to enqueue tasks
    template <class F>
    void enqueue(F&& f);

    // Non-template methods
    void restart(size_t numThreads);
    void stopPool();

private:
    void start(size_t numThreads);

    std::vector<std::thread> workers;
    std::queue<std::function<void()>> taskQueue;
    std::mutex queueMutex;
    std::condition_variable condition;
    bool stop;
};

// Template method definition directly in the header
template <class F>
void ThreadPool::enqueue(F&& f) {
    {
        std::unique_lock<std::mutex> lock(queueMutex);
        if (stop) {
            throw std::runtime_error("ThreadPool is stopped, cannot enqueue tasks");
        }
        taskQueue.push(std::forward<F>(f)); // Forward the task to the queue
    }
    condition.notify_one();
}

#endif