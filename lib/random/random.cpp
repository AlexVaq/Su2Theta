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
	double	max     () { return (*myRNG).Max(); }
	double	min     () { return (*myRNG).Min(); }

//	template<class T>
//	T	genVRand() {
//		return	T(0.);
//	}

	template<>
	float	genVRand<float>() {
		return	(*myRNG)();
	}

	template<>
	double	genVRand<double>() {
		return	(*myRNG)();
	}

	template<>
	float	vGenRand<float>() {
		return	(*myRNG)();
	}

	template<>
	double	vGenRand<double>() {
		return	(*myRNG)();
	}

	template<>
	Simd::Simd_f	genVRand<Simd::Simd_f>() {
		float in[Simd::Simd_f::sWide];

		#pragma unroll
		for (int i=0; i<Simd::Simd_f::sWide; i++)
			in[i] = genRand();

		return	Simd::Simd_f(in);
	}

	template<>
	Simd::Simd_d	genVRand<Simd::Simd_d>() {
		double in[Simd::Simd_d::sWide];

		#pragma unroll
		for (int i=0; i<Simd::Simd_d::sWide; i++)
			in[i] = genRand();

		return	Simd::Simd_d(in);
	}
}
