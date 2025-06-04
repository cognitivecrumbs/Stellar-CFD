#ifndef WORKER_CLASS_H
#define WORKER_CLASS_H

#include "ThreadPool.h"

class WorkerClass {
public:
    explicit WorkerClass(ThreadPool& threadPool);
    void executeTasks();

private:
    ThreadPool& pool;
};

#endif