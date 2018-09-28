#include "simd/simd.h"
#include "random/random.h"

namespace	Simd {

	Simd_f  sqrt    (const Simd_f &x) {
                return	opCode(sqrt_ps, x.data);
        }

        Simd_f  cos     (const Simd_f &x) {
                return	opCode(cos_ps, x.data);
        }

        Simd_f  sin     (const Simd_f &x) {
                return	opCode(sin_ps, x.data);
        }

        Simd_f  log     (const Simd_f &x) {
                return	opCode(log_ps, x.data);
        }

        Simd_f  exp	(const Simd_f &x) {
                return	opCode(exp_ps, x.data);
        }

        Simd_f  abs	(const Simd_f &x) {
		#ifdef	__AVX512F__
                	return	opCode(abs_ps, x.data);
		#else
                	return	opCode(andnot_ps, opCode(castsi256_ps, iSgnAbsf), x.data);
		#endif
        }

	inline	void	Simd_f::SetRandom () {
		#ifdef	__AVX512F__
			(*this) = Simd_f(Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand(),
					 Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand(),
					 Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand(),
					 Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand());
		#else
			(*this) = Simd_f(Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand(),
					 Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand());
		#endif
	}

	Simd_d  sqrt	(const Simd_d &x) {
                return	opCode(sqrt_pd, x.data);
        }

        Simd_d  cos	(const Simd_d &x) {
                return	opCode(cos_pd, x.data);
        }

        Simd_d  sin	(const Simd_d &x) {
                return	opCode(sin_pd, x.data);
        }

        Simd_d  log	(const Simd_d &x) {
                return	opCode(log_pd, x.data);
        }

        Simd_d  exp	(const Simd_d &x) {
                return	opCode(exp_pd, x.data);
        }

        Simd_d  abs	(const Simd_d &x) {
		#ifdef	__AVX512F__
                	return	opCode(abs_pd, x.data);
		#else
                	return	opCode(andnot_pd, opCode(castsi256_pd, iSgnAbsd), x.data);
		#endif
        }

	inline	void	Simd_d::SetRandom () {
		#ifdef	__AVX512F__
			(*this) = Simd_d(Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand(),
					 Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand());
		#else
			(*this) = Simd_d(Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand());
		#endif
	}
};
