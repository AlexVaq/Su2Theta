#ifndef	__COMMS
	#define	__COMMS

	#include <climits>
	#include "enumFields.h"

	#ifndef __USE_POSIX
		#define HOST_NAME_MAX   128
	#endif

	namespace	Su2Comms {
		extern int  rank;
		extern int  ranksPerNode;
		extern int  nThreads;
		extern int  idxAcc;
		extern int  size;
		extern char hostname[HOST_NAME_MAX];

		extern int  maxThreadsPerBlock;
		extern int  maxThreadsPerDim[3];
		extern int  maxGridSize[3];

		extern size_t gpuMem;

		int	initComms (int argc, char *argv[], int size, Su2Enum::Device dev, Su2Enum::LogMpi logMpi, Su2Enum::VerbosityLevel verb);
		void	endComms();
		void	commSync();
	}
#endif
