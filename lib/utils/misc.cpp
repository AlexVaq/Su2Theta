#include "utils/logger.h"
#include "utils/memAlloc.h"
#include "utils/profiler.h"
#include "random/random.h"
//#include "utils/parse.h"

using namespace Su2Prof;
using namespace Su2Comms;

void	initSu2 (int argc, char *argv[]) {
//	parseArgs(argc, argv);

	if (Su2Comms::initComms(argc, argv, 1, DeviceCpu, ZeroRank, VerbHigh) == -1)
	{
		LogOut ("Error initializing devices and Mpi\n");
        	exit(1);
	}

//	if (Su2Comms::rank == 0)
//		createOutput();

//	LogMsg (VERB_NORMAL, "Output folder set to %s", outDir);
//	LogMsg (VERB_NORMAL, "FFTW wisdom folder set to %s", wisDir);

	Su2Prof::initProfilers();
	Su2Rand::initRandom();

	return;
}

void	endSu2() {
	printMemStats();
	printProfStats();
	Su2Comms::endComms();

	exit(0);

	return;
}
