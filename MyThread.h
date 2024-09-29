#ifndef SUBTHREAD_H
#define SUBTHREAD_H

#include "threads.h"

using namespace std;

//
class MyThread : public Thread
{
public:
	// This method will be executed by the Thread::exec method,
	// which is executed in the thread of execution
	virtual void run();
};

#endif
