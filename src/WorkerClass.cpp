#include <WorkerClass.h>
#include <iostream>
#include <thread>
#include <chrono>

WorkerClass::WorkerClass(ThreadPool& threadPool) : pool(threadPool) {}

void WorkerClass::executeTasks() {
    pool.enqueue([] {
        std::cout << "Task 1: Running from WorkerClass!" << std::endl;
    });

    pool.enqueue([] {
        std::cout << "Task 2: Doing some work..." << std::endl;
        std::this_thread::sleep_for(std::chrono::seconds(1));
        std::cout << "Task 2: Work done!" << std::endl;
    });
}