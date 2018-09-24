#include "simd/simd.h"

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
    /* 
     	Simd_f	log2	(const Simd_f &x) {
                return	opCode(set_ps, std::log(x.data[7]), std::log(x.data[6]), std::log(x.data[5]), std::log(x.data[4]), std::log(x.data[3]), std::log(x.data[2]), std::log(x.data[1]), std::log(x.data[0]));
        }

        Simd_f	exp2	(const Simd_f &x) {
                return	opCode(set_ps, std::exp(x.data[7]), std::exp(x.data[6]), std::exp(x.data[5]), std::exp(x.data[4]), std::exp(x.data[3]), std::exp(x.data[2]), std::exp(x.data[1]), std::exp(x.data[0]));
        }
*/

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
/*
        Simd_d  log2	(const Simd_d &x) {
                return	opCode(set_pd, std::log(x.data[3]), std::log(x.data[2]), std::log(x.data[1]), std::log(x.data[0]));
        }

        Simd_d  exp2	(const Simd_d &x) {
                return	opCode(set_pd, std::exp(x.data[3]), std::exp(x.data[2]), std::exp(x.data[1]), std::exp(x.data[0]));
        }
*/
        Simd_d  abs	(const Simd_d &x) {
		#ifdef	__AVX512F__
                	return	opCode(abs_pd, x.data);
		#else
                	return	opCode(andnot_pd, opCode(castsi256_pd, iSgnAbsd), x.data);
		#endif
        }
};
