// #include "ThreadPool.h"
// #include "WorkerClass.h"

#include <ThreadPool.h>
#include <WorkerClass.h>
#include <thread>
#include <chrono>

int main() {
    ThreadPool pool(4);
    WorkerClass worker(pool);

    worker.executeTasks();

    // Allow time for tasks to complete
    std::this_thread::sleep_for(std::chrono::seconds(2));
    return 0;
}