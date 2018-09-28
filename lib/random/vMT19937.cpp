#include <vector>
#include <random>
#include <memory>

#include "enumFields.h"
#include "simd/simd.h"
#include "random/random.h"
#include "utils/logger.h"

#include <omp.h>

using namespace	Su2Enum;
using namespace	Simd;

namespace	Su2Rand {

	/*	We need two different generators for different vector lengths	*/
	std::shared_ptr<vRandomGen<Simd_f>> vRNG_f;
	std::shared_ptr<vRandomGen<Simd_d>> vRNG_d;

	void	vInitRandom () {
		static bool	initRNG = false;

		if (initRNG == false) {
			vRNG_f = std::make_shared<vRandomGen<Simd_f>>();
			vRNG_d = std::make_shared<vRandomGen<Simd_d>>();
		} else {
			LogMsg(VerbHigh, "Vectorized random number generator already initialized");
		}
	}

	template<>
	Simd_f	vGenRand<Simd_f>() {
		return	(*vRNG_f)();
	}

	template<>
	Simd_d	vGenRand<Simd_d>() {
		return	(*vRNG_d)();
	}

	template<>
	void	RNG<Simd_f>::initRNG(uint seed) {

		uint	*iState  = static_cast<uint  *>  (static_cast<void*>(state));

		iState[0] = seed;
		iState[1] = mtInit * (iState[0]^(iState[0] >> 30)) + 1;
		iState[2] = mtInit * (iState[1]^(iState[1] >> 30)) + 2;
		iState[3] = mtInit * (iState[2]^(iState[2] >> 30)) + 3;

		for (int i=4; i<Simd_f::sWide; i+=4) {
			iState[i+0] = ((104729 + iState[i-4])*691) % 179426549;
			iState[i+1] = mtInit * (iState[i+0]^(iState[i+0] >> 30)) + 1;
			iState[i+2] = mtInit * (iState[i+1]^(iState[i+1] >> 30)) + 2;
			iState[i+3] = mtInit * (iState[i+2]^(iState[i+2] >> 30)) + 3;
		}

		for (int i=Simd_f::sWide; i<(mtSize32 + 1)*Simd_f::sWide; i+=Simd_f::sWide)
			for (int j=0; j<Simd_f::sWide; j+=4) {
				iState[i+j+0] = mtInit * (iState[i+j+3-Simd_f::sWide]^(iState[i+j+3-Simd_f::sWide] >> 30)) + (i<<2)/Simd_f::sWide;
				iState[i+j+1] = mtInit * (iState[i+j+0]^(iState[i+j+0] >> 30)) + (i<<2)/Simd_f::sWide + 1;
				iState[i+j+2] = mtInit * (iState[i+j+1]^(iState[i+j+1] >> 30)) + (i<<2)/Simd_f::sWide + 2;
				iState[i+j+3] = mtInit * (iState[i+j+2]^(iState[i+j+2] >> 30)) + (i<<2)/Simd_f::sWide + 3;
			}

//		for (int i=0; i<mtSize32*Simd_f::sWide; i++)
//			iState[i] = (iState[i] & mtLowMask32) | mtHighConst32;

		/*	Period certification	*/
		for (int i=0; i<Simd_f::sWide; i+=4) {
			auto inner = ((iState[i] & mtPcf1) ^ (iState[i+3] & mtPcf2));

			for (int i=16; i>0; i>>=1)
				inner ^= inner >> i;

			inner &= 1;

			if (inner != 1)
				iState[i] ^= 1;
		}

		idx = mtSize32;
	}
                           
	template<>         
	void	RNG<Simd_d>::initRNG(uint seed) {

		uint	*iState  = static_cast<uint  *>  (static_cast<void*>(state));
		uint64	*iSta64  = static_cast<uint64*>  (static_cast<void*>(state));

		iState[0] = seed;
		iState[1] = mtInit * (iState[0]^(iState[0] >> 30)) + 1;

		for (int i=2; i<Simd_d::sWide; i+=2) {
			iState[i+0] = ((104729 + iState[i-2])*691) % 179426549;
			iState[i+1] = mtInit * (iState[i]^(iState[i] >> 30)) + 1;
		}

		for (int i=Simd_d::sWide; i<(mtSize64 + 1)*Simd_d::sWide; i+=Simd_d::sWide)
			for (int j=0; j<Simd_d::sWide; j+=2) {
				//iState[i+j+0] = mtInit * (iState[i+j+1-Simd_d::sWide]^(iState[i+j+1-Simd_d::sWide] >> 30)) + (i<<1)/Simd_d::sWide + (j%2);
				iState[i+j+0] = mtInit * (iState[i+j+3-Simd_d::sWide]^(iState[i+j+3-Simd_d::sWide] >> 30)) + (i<<2)/Simd_d::sWide;
				iState[i+j+1] = mtInit * (iState[i+j+0]^(iState[i+j+0] >> 30)) + (i<<1)/Simd_d::sWide + 1;
			}

		for (int i=0; i<mtSize64*Simd_d::sWide; i++)
			iSta64[i] = (iSta64[i] & mtLowMask64) | mtHighConst64;

		/*	Period certification	*/
		for (int i=0; i<Simd_d::sWide; i+=2) {
			auto inner = ((iSta64[i] ^ mtFix1) & mtPcd1) ^ ((iSta64[i+1] ^ mtFix2) & mtPcd2);

			for (int i=32; i>0; i>>=1)
				inner ^= inner >> i;

			inner &= 1;

			if (inner != 1)
				iSta64[i+1] ^= 1;
		}

		idx = mtSize64;
	}

	template<>
	void	RNG<Simd_f>::genNext() {

		Simd_f c = state[mtSize32-2];
		Simd_f d = state[mtSize32-1];

		for (int i=0; i<mtSize32; i++) {
			auto a = state[i];
			auto x = a.glShift(1);
			auto y = (state[(i+mtPos32)%mtSize32] >> mtSr32) & mtMsk32;
			auto z = c.grShift(1);
			auto v = d << mtSl32;

			c = d;
			d = (((z ^ a) ^ v) ^ x) ^ y;
			state[i] = d;
		}
	}
                           
	template<>
	void	RNG<Simd_d>::genNext() {

		Simd_d u = state[mtSize64];

		for (int i=0; i<mtSize64; i++) {
			auto x = state[i];
			auto z = x << mtSl64;
			auto y = u.rPermute();

			z ^= state[(i+mtPos64)%mtSize64];
			y ^= z;

			y >>= mtSr64;
			state[i] = y ^ x;
			u = y & mtMsk64;
		}

		state[mtSize64] = u;
	}

	template<>
	Simd_f	RNG<Simd_f>::operator()() {
		if (idx >= mtSize32) {
			// Regenerate a bunch of random numbers
			genNext();
			idx = 0;
		}

		return	((state[idx++] & mtLow32) | mtHigh32) - Simd_f(1.0);
	}

	template<>
	Simd_d	RNG<Simd_d>::operator()() {
		if (idx >= mtSize64) {
			// Regenerate a bunch of random numbers
			genNext();
			idx = 0;
		}

		return	state[idx++] - Simd_d(1.0);
	}

/*
	template<class T>
	T	RNG<T>::operator()() {
		if (idx >= mtSize) {
			// Regenerate a bunch of random numbers
			genNext();
			idx = 0;
		}

		return	state[idx++] - (T)(1.0);
	}
*/
}
