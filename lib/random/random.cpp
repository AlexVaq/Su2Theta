#include <vector>
#include <random>
#include <memory>

#include "enumFields.h"
#include "random/random.h"
#include "utils/logger.h"

#include <omp.h>

using namespace	Su2Enum;

namespace	Su2Rand {

	std::shared_ptr<RandomGen<double>> myRNG;

	void	initRandom () {
		static bool	initRNG = false;

		if (initRNG == false) {
			myRNG = std::make_shared<RandomGen<double>>();
		} else {
			LogMsg(VerbHigh, "Random number generator already initialized");
		}
	}

	double	genRand () { return (*myRNG)(); }
}
