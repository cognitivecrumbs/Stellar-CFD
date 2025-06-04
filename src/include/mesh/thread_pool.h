#ifndef THREADPOOL_H
#define THREADPOOL_H
#include <iostream>
#include <vector>
#include <thread>
#include <functional>
#include <queue>
#include <stack>
#include <mutex>
#include <condition_variable>
#include <atomic>

// yanked code, not sure how well it works

// class ThreadPool {
// public:
//     ThreadPool(size_t numThreads) : stop(false) {
//         for (size_t i = 0; i < numThreads; ++i) {
//             workers.emplace_back([this] {
//                 while (true) {
//                     std::function<void()> task;
//                     {
//                         std::unique_lock<std::mutex> lock(this->queueMutex);
//                         this->condition.wait(lock, [this] { return this->stop || !this->taskQueue.empty(); });
//                         if (this->stop && this->taskQueue.empty()) {
//                             return;
//                         }
//                         task = std::move(this->taskQueue.front());
//                         this->taskQueue.pop();
//                     }
//                     task();
//                 }
//             });
//         }
//     }

//     template <class F>
//     void enqueue(F&& f) {
//         {
//             std::unique_lock<std::mutex> lock(queueMutex);
//             if (stop) {
//                 throw std::runtime_error("ThreadPool is stopped, cannot enqueue tasks");
//             }
//             taskQueue.push(std::forward<F>(f));
//         }
//         condition.notify_one();
//     }

//     void stopPool() {
//         {
//             std::unique_lock<std::mutex> lock(queueMutex);
//             stop = true;
//         }
//         condition.notify_all();
//         for (std::thread &worker : workers) {
//             worker.join();
//         }
//     }

//     ~ThreadPool() {
//         if (!stop) {
//             stopPool();
//         }
//     }

// private:
//     std::vector<std::thread> workers;
//     std::queue<std::function<void()>> taskQueue;
//     std::mutex queueMutex;
//     std::condition_variable condition;
//     bool stop;
// };

// #include <vector>
// #include <queue>
// #include <thread>
// #include <mutex>
// #include <condition_variable>
// #include <future>
// #include <functional>
// #include <stdexcept>

// class ThreadPool {
// public:
//     ThreadPool(size_t numThreads) : stop(false) {
//         start(numThreads);
//     }

//     // Method to restart the thread pool
//     void restart(size_t numThreads) {
//         stopPool();  // Stop any existing threads
//         start(numThreads);  // Start the thread pool again with a new number of threads
//     }

//     template <class F, class... Args>
//     auto submit(F&& f, Args&&... args) -> std::future<decltype(f(args...))> {
//         using ReturnType = decltype(f(args...));

//         auto task = std::make_shared<std::packaged_task<ReturnType()>>(
//             std::bind(std::forward<F>(f), std::forward<Args>(args)...)
//         );

//         std::future<ReturnType> result = task->get_future();
//         {
//             std::unique_lock<std::mutex> lock(queueMutex);
//             if (stop) {
//                 throw std::runtime_error("ThreadPool is stopped, cannot submit tasks");
//             }
//             taskQueue.emplace([task]() { (*task)(); });
//         }
//         condition.notify_one();
//         return result;
//     }

//     void stopPool() {
//         {
//             std::unique_lock<std::mutex> lock(queueMutex);
//             stop = true;
//         }
//         condition.notify_all();
//         for (std::thread &worker : workers) {
//             worker.join();
//         }
//         workers.clear();  // Clear the worker threads after stopping
//     }

//     ~ThreadPool() {
//         if (!stop) {
//             stopPool();
//         }
//     }

// private:
//     void start(size_t numThreads) {
//         stop = false;
//         workers.reserve(numThreads);
//         for (size_t i = 0; i < numThreads; ++i) {
//             workers.emplace_back([this] {
//                 while (true) {
//                     std::function<void()> task;
//                     {
//                         std::unique_lock<std::mutex> lock(this->queueMutex);
//                         this->condition.wait(lock, [this] { return this->stop || !this->taskQueue.empty(); });
//                         if (this->stop && this->taskQueue.empty()) {
//                             return;
//                         }
//                         task = std::move(this->taskQueue.front());
//                         this->taskQueue.pop();
//                     }
//                     task();
//                 }
//             });
//         }
//     }

//     std::vector<std::thread> workers;
//     std::queue<std::function<void()>> taskQueue;
//     std::mutex queueMutex;
//     std::condition_variable condition;
//     bool stop;
// };

// --------------------------------

class ThreadPool {
public:
    ThreadPool(size_t numThreads) : stop(false) {
        start(numThreads);
    }

    // Method to restart the thread pool
    void restart(size_t numThreads) {
        stopPool();  // Stop any existing threads
        start(numThreads);  // Start the thread pool again with a new number of threads
    }

    template <class F>
    void enqueue(F&& f) {
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            if (stop) {
                throw std::runtime_error("ThreadPool is stopped, cannot enqueue tasks");
            }
            taskQueue.push(std::forward<F>(f));
        }
        condition.notify_one();
    }

    void stopPool() {
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            stop = true;
        }
        condition.notify_all();
        for (std::thread &worker : workers) {
            worker.join();
        }
        workers.clear();  // Clear the worker threads after stopping
    }
    
    // Running stop pool: allow submitting new tasks but waits until all tasks are completed before stopping
    void runningStopPool() {
        
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            // stop = true;
        }
        condition.notify_all();  // Wake up workers to check the stop condition

        // Wait until all tasks are finished before truly stopping the pool
        while (activeTasks > 0) {
            std::this_thread::yield();  // Yield the thread to let other threads finish their tasks
        }

        // Join all worker threads
        for (std::thread& worker : workers) {
            worker.join();
        }
        workers.clear();  // Clear the worker threads after stopping
    }

    ~ThreadPool() {
        if (!stop) {
            stopPool();
        }
    }

private:
    void start(size_t numThreads) {
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

    // void start(size_t numThreads) {
    //     workers.reserve(numThreads);
    //     for (size_t i = 0; i < numThreads; ++i) {
    //         workers.emplace_back([this] {
    //             while (true) {
    //                 std::function<void()> task;
    //                 {
    //                     std::unique_lock<std::mutex> lock(this->queueMutex);
    //                     this->condition.wait(lock, [this] {
    //                         return this->stop || !this->taskQueue.empty();
    //                     });

    //                     // Exit if stop is true and no tasks are left
    //                     if (this->stop && this->taskQueue.empty()) {
    //                         return;
    //                     }

    //                     task = std::move(this->taskQueue.top());
    //                     this->taskQueue.pop();
    //                 }

    //                 // Increment active task count
    //                 {
    //                     std::unique_lock<std::mutex> lock(this->queueMutex);
    //                     ++activeTasks;
    //                 }

    //                 task();  // Execute the task

    //                 // Decrement active task count after the task is done
    //                 {
    //                     std::unique_lock<std::mutex> lock(this->queueMutex);
    //                     --activeTasks;
    //                 }
    //             }
    //         });
    //     }
    // }

    // void start(size_t numThreads) {
    //     workers.reserve(numThreads);
    //     for (size_t i = 0; i < numThreads; ++i) {
    //         workers.emplace_back([this] {
    //             while (true) {
    //                 std::function<void()> task;
    //                 {
    //                     std::unique_lock<std::mutex> lock(this->queueMutex);
                        
    //                     // Pause if there are no tasks, and wait for a new one
    //                     this->condition.wait(lock, [this] {
    //                         return this->stop || !this->taskQueue.empty();
    //                     });

    //                     // Exit if stop is true and no tasks are left
    //                     if (this->stop && this->taskQueue.empty()) {
    //                         return;
    //                     }

    //                     // Get the latest task (LIFO from stack)
    //                     task = std::move(this->taskQueue.top());
    //                     this->taskQueue.pop();
    //                 }

    //                 // Increment active task count
    //                 {
    //                     std::unique_lock<std::mutex> lock(this->queueMutex);
    //                     ++activeTasks;
    //                 }

    //                 // Execute the task
    //                 task();  

    //                 // Decrement active task count after the task is done
    //                 {
    //                     std::unique_lock<std::mutex> lock(this->queueMutex);
    //                     --activeTasks;
    //                 }
    //             }
    //         });
    //     }
    // }

    std::vector<std::thread> workers;
    std::queue<std::function<void()>> taskQueue;
    // std::stack<std::function<void()>> taskQueue;
    std::mutex queueMutex;
    std::condition_variable condition;
    std::atomic<int> activeTasks;

    bool stop;
};

// #include <iostream>
// #include <vector>
// #include <queue>
// #include <thread>
// #include <mutex>
// #include <condition_variable>
// #include <functional>
// #include <atomic>
// #include <stdexcept>

// class ThreadPool {
// public:
//     ThreadPool(size_t numThreads) : stop(false), activeTasks(0) {
//         start(numThreads);
//     }

//     ~ThreadPool() {
//         stopPool();
//     }

//     // Method to restart the thread pool with a new number of threads
//     void restart(size_t numThreads) {
//         stopPool();  // Stop existing threads
//         start(numThreads);  // Start with new threads
//     }

//     // Method to enqueue a task
//     template <class F>
//     void enqueue(F&& f) {
//         {
//             std::unique_lock<std::mutex> lock(queueMutex);
//             if (stop) {
//                 throw std::runtime_error("ThreadPool is stopped, cannot enqueue tasks");
//             }
//             taskQueue.push(std::forward<F>(f));
//         }
//         condition.notify_one();  // Notify a worker to process the task
//     }

//     // Stop the pool after all tasks are finished
//     void stopPool() {
//         {
//             std::unique_lock<std::mutex> lock(queueMutex);
//             stop = true;
//         }
//         condition.notify_all();  // Wake up all workers to allow them to stop

//         // Wait until all tasks are finished
//         while (activeTasks > 0) {
//             std::this_thread::yield();  // Yield to allow other threads to finish their tasks
//         }

//         // Join all threads
//         for (std::thread& worker : workers) {
//             worker.join();
//         }
//         workers.clear();  // Clear the worker threads after stopping
//     }

// private:
//     void start(size_t numThreads) {
//         workers.reserve(numThreads);
//         for (size_t i = 0; i < numThreads; ++i) {
//             workers.emplace_back([this] {
//                 while (true) {
//                     std::function<void()> task;
//                     {
//                         std::unique_lock<std::mutex> lock(this->queueMutex);
//                         this->condition.wait(lock, [this] {
//                             return this->stop || !this->taskQueue.empty();
//                         });

//                         if (this->stop && this->taskQueue.empty()) {
//                             return;  // Exit if pool is stopped and no tasks are left
//                         }

//                         task = std::move(this->taskQueue.front());
//                         this->taskQueue.pop();
//                     }

//                     // Increment active task count
//                     {
//                         std::unique_lock<std::mutex> lock(this->queueMutex);
//                         ++activeTasks;
//                     }

//                     task();  // Execute the task

//                     // Decrement active task count after the task is done
//                     {
//                         std::unique_lock<std::mutex> lock(this->queueMutex);
//                         --activeTasks;
//                     }
//                 }
//             });
//         }
//     }

//     std::vector<std::thread> workers;
//     std::queue<std::function<void()>> taskQueue;
//     std::mutex queueMutex;
//     std::condition_variable condition;
//     bool stop;

//     // Atomic counter to track active tasks
//     std::atomic<int> activeTasks;
// };

#endif // THREADPOOL_H