#ifndef	__COMMS
	#define	__COMMS

	#include <cstddef>
	#include <climits>

	#ifndef __USE_POSIX
		#define HOST_NAME_MAX   128
	#endif

	using namespace	std;

	namespace	Su2Comms {
		int  rank          =  0;
		int  ranksPerNode  =  1;
		int  nThreads      =  1;
		int  idxAcc        = -1;
		int  size          =  1;
		char hostname[HOST_NAME_MAX];

		int  maxThreadsPerBlock  = 0;
		int  maxThreadsPerDim[3] = { 0, 0, 0 };
		int  maxGridSize[3]      = { 0, 0, 0 };

		size_t gpuMem            = 0;
	}
#endif
