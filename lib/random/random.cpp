#include <memory>

#include "enumFields.h"
#include "random/random.h"
#include "comms/comms.h"
#include "utils/logger.h"

#include <omp.h>

using namespace	Su2Enum;

namespace	Su2Rand {

	template<typename Float>
		RandomGen<Float>::RandomGen() {

			decltype(uni.param()) newRange (0.0, 1.0);

			sd.resize  (Su2Comms::nThreads);
			mt64.resize(Su2Comms::nThreads);
			uni.param  (newRange);

			std::random_device seed;			// Totally random seed coming from memory garbage

			for (int i=0; i<Su2Comms::nThreads; i++)
				sd[i] = seed()*(1 + Su2Comms::rank);

			#pragma omp parallel default(shared)
			{
				int nThread = omp_get_thread_num();
				mt64[nThread].seed(sd[nThread]);	// Mersenne-Twister 64 bits, independent per thread
			}
		}

	template<typename Float>
		RandomGen<Float>::RandomGen(int seed) {

			decltype(uni.param()) newRange (0.0, 1.0);

			sd.resize  (Su2Comms::nThreads);
			mt64.resize(Su2Comms::nThreads);
			uni.param  (newRange);

			for (int i=0; i<Su2Comms::nThreads; i++)
				sd[i] = seed*(1 + Su2Comms::rank);


			#pragma omp parallel default(shared)
			{
				int nThread = omp_get_thread_num();
				mt64[nThread].seed(sd[nThread]);	// Mersenne-Twister 64 bits, independent per thread
			}
		}


	template<typename Float>
	Float	RandomGen<Float>::rand() {
		int nThread = omp_get_thread_num();
		return	uni(mt64[nThread]);
	}

	std::shared_ptr<RandomGen<double>> myRNG;

	void	initRandom () {
		static bool	initRNG = false;

		if (initRNG == false) {
			myRNG = std::make_shared<RandomGen<double>>();
		} else {
			LogMsg(VerbHigh, "Random number generator already initialized");
		}
	}
}
