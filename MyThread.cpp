//==================================================================
// File		: mythread.h
// Author	: rsagalyn
// Date		: Aug 25, 2013
// Description	: Subclass of threads.h, implements run()
//==================================================================

#include "MyThread.h"

void MyThread::run()
{

	if (algorithm == TSP::Algorithm::DP)
		mytsp->openTSP_DP(my_id);
	else
		mytsp->openTSP_EMST(my_id);

	// sleep for a bit

	// if (DEBUG) cout << "thread " << setw(4) << left << my_id << setw(8) << left
	// 		<< " result: " << setw(5) << left << result << endl;

	pthread_exit(NULL);
}
