#ifndef	__RANDOMGEN
	#define	__RANDOMGEN
	#include <vector>
	#include <random>
	#include <memory>

	#include "enumFields.h"
	#include "comms/comms.h"
	#include "utils/logger.h"

	#include <omp.h>

	using namespace	Su2Enum;

	namespace	Su2Rand {

		void	initRandom ();

		template<typename Float>
		class	RandomGen {
			private:

			std::vector<int>			sd;
			std::vector<std::mt19937_64>		mt64;
			std::uniform_real_distribution<Float>	uni;

			public:

				RandomGen() {

				decltype(uni.param()) newRange (0.0, 1.0);

				sd.resize  (Su2Comms::nThreads);
				mt64.resize(Su2Comms::nThreads);
				uni.param  (newRange);

				std::random_device seed;			// Totally random seed coming from memory garbage

				for (int i=0; i<Su2Comms::nThreads; i++)
					sd[i] = i;//seed()*(1 + Su2Comms::rank);

				#pragma omp parallel default(shared)
				{
					int nThread = omp_get_thread_num();
					mt64[nThread].seed(sd[nThread]);	// Mersenne-Twister 64 bits, independent per thread
				}
			}

				RandomGen(int seed) {

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

			Float	operator()() {
				int nThread = omp_get_thread_num();
				return	uni(mt64[nThread]);
			}

			double	Min() {
				return	uni.min();
			}

			double	Max() {
				return	uni.max();
			}
		};

		double	genRand	();
		double	max	();
		double	min	();

		template<class T>
		T	genVRand();
	}

#include <random/vMT19937.h>
#endif
