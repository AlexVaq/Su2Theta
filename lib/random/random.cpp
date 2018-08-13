#include <vector>
#include <random>
#include <memory>

#include "enumFields.h"
#include "random/random.h"
#include "simd/simd.h"
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

//	template<class T>
//	T	genVRand() {
//		return	T(0.);
//	}

	template<>
	Simd::Simd_f	genVRand<Simd::Simd_f>() {
		float in[Simd::Simd_f::nData];

		#pragma unroll
		for (int i=0; i<Simd::Simd_f::nData; i++)
			in[i] = genRand();

		return	Simd::Simd_f(in);
	}
}
