#ifndef	__LOGEXPINTRINSICS
#define	__LOGEXPINTRINSICS

#include<cmath>

#define opCode_P(x,y,...) x ## _ ## y (__VA_ARGS__)
#define opCode_N(x,y,...) opCode_P(x, y, __VA_ARGS__)
#define opCode(x,...) opCode_N(_PREFIX_, x, __VA_ARGS__)

#include <immintrin.h>

#ifdef	__AVX512F__
	#define _MData_ __m512d
	#define	_MInt_  __m512i
	#define	_MHnt_  __m256i
#elif   defined(__AVX__)
	#define _MData_ __m256d
	#define	_MInt_  __m256i
	#define	_MHnt_  __m128i
#else
	#define _MData_ __m128d
	#define	_MInt_  __m128i
#endif

#if	defined(__AVX512F__)
	#define	_PREFIX_ _mm512
	#define	_PREFXL_ _mm256
	#define opCodl(x,...) opCode_N(_PREFXL_, x, __VA_ARGS__)
#else
	#if not defined(__AVX__) and not defined(__AVX2__)
		#define	_PREFIX_ _mm
	#else
		#define	_PREFIX_ _mm256
		#define	_PREFXL_ _mm
		#define opCodl(x,...) opCode_N(_PREFXL_, x, __VA_ARGS__)
	#endif
#endif

#define	M_PI2	(M_PI *M_PI)
#define	M_PI4	(M_PI2*M_PI2)
#define	M_PI6	(M_PI4*M_PI2)

#ifndef	__INTEL_COMPILER
#if 0
inline _MData_	opCode(exp2_pd, _MData_ &x) {
	_MInt_	N, n1, n2, M;
	_MData_ n1f, n2f, R1, R2, R, Q, Sl, St, S;
/*
	1.  Filtra Nan
	2.  Filtra +inf -> +inf
	3.  Filtra -inf -> 0
        4.  Filtra Threshold_1 -> +inf / 0
	5.  Filtra Threshold_2 -> 1+x
*/

	N   = opCode(cvtps_epi32, opCode(round_ps, opCode(mul_ps, x, vIvL_f), _MM_FNEAREST_TO_NEAREST_INT | _MM_FROUND_NO_EXC));
	n2  = opCode(and_si256, N, m32Mask);	// Module 32
	n1  = opCode(sub_epi32, N, n2);
	n1f = opCode(cvtepi32_ps, n1)
        n2f = opCode(cvtepi32_ps, n2)
	R1  = opCode(sub_ps,
		opCode(sub_ps, x, opCode(mul_ps, n1f, vL1f)),
		opCode(mul_ps, n2f, vL1f));
	// 9.  R1 = (x - n1*L1) - n2*L1 / R1 = (x - N*L1)
	R2  = opCode(mul_ps, N, vL2f);	// - Sign in vL2f
	M   = opCode(srli_epi32, n1, 5);
	R   = opCode(add_ps, R1, R2);
	Q   = opCode(mul_ps,
		opCode(mul_ps, R, R),
		opCode(add_ps, vA1f, opCode(mul_ps, R, vA2f)));
	Q   = opCode(add_ps, R1, opCode(add_ps, R2, Q));
	Sl  = opCode(set1r_ps, sLead_f[n2[7]], sLead_f[n2[6]], sLead_f[n2[5]], sLead_f[n2[4]], sLead_f[n2[3]], sLead_f[n2[2]], sLead_f[n2[1]], sLead_f[n2[0]]);
	St  = opCode(set1r_ps, sTrail_f[n2[7]], sTrail_f[n2[6]], sTrail_f[n2[5]], sTrail_f[n2[4]], sTrail_f[n2[3]], sTrail_f[n2[2]], sTrail_f[n2[1]], sTrail_f[n2[0]]);
	S   = opCode(add_ps, Sl, St);
	return	opCode(mul_ps,
			opCode(cvtepi32_ps, opCode(sllv_epi32, two, M)),
			opCode(add_ps, Sl, opCode(add_ps, St, opCode(mul_ps, S, Q))));
}
#endif
#ifdef	__AVX512F__
inline _MData_	opCode(exp_pd, _MData_ &x)
{
	return	opCode(set_pd, std::exp(x[7]), std::exp(x[6]), std::exp(x[5]), std::exp(x[4]), std::exp(x[3]), std::exp(x[2]), std::exp(x[1]), std::exp(x[0]));
}
#elif   defined(__AVX__)
inline _MData_	opCode(exp_pd, _MData_ &x)
{
	return	opCode(set_pd, std::exp(x[3]), std::exp(x[2]), std::exp(x[1]), std::exp(x[0]));
}
#else
inline _MData_	opCode(exp_pd, _MData_ &x)
{
	return	opCode(set_pd, std::exp(x[1]), std::exp(x[0]));
}
#endif

#endif

#undef	_MData_

#if	defined(__AVX512F__)
	#define	_MData_ __m512
#elif	defined(__AVX__)
	#define	_MData_ __m256
#else
	#define	_MData_ __m128
#endif 

#ifndef	__INTEL_COMPILER

inline _MData_	opCode(exp2_ps, _MData_ &x) {
	_MInt_	N, n1, n2, M;
	_MData_ n1f, n2f, R1, R2, R, Q, Sl, St, S;

/*
	1.  Filtra Nan
	2.  Filtra +inf -> +inf
	3.  Filtra -inf -> 0
        4.  Filtra Threshold_1 -> +inf / 0
	5.  Filtra Threshold_2 -> 1+x	<-- No hace falta, porque vamos a calcular la chunga igual...
*/

	N   = opCode(cvtps_epi32, opCode(round_ps, opCode(mul_ps, x, vIvL_f), _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC));
	n2  = opCode(and_si256, N, m32Mask);	// Module 32
	n1  = opCode(sub_epi32, N, n2);
	n1f = opCode(cvtepi32_ps, n1);
        n2f = opCode(cvtepi32_ps, n2);
	R1  = opCode(sub_ps,
		opCode(sub_ps, x, opCode(mul_ps, n1f, vL1f)),
		opCode(mul_ps, n2f, vL1f));
	// 9.  R1 = (x - n1*L1) - n2*L1 / R1 = (x - N*L1)
	R2  = opCode(mul_ps, opCode(cvtepi32_ps, N), vL2f);	// - Sign in vL2f
	M   = opCode(srli_epi32, n1, 5);
	R   = opCode(add_ps, R1, R2);
	Q   = opCode(mul_ps,
		opCode(mul_ps, R, R),
		opCode(add_ps, vA1f, opCode(mul_ps, R, vA2f)));
	Q   = opCode(add_ps, R1, opCode(add_ps, R2, Q));
	int vals[8];
	opCode(store_si256, static_cast<__m256i*>(static_cast<void*>(vals)), n2);
	Sl  = opCode(set_ps, sLead_f[vals[7]], sLead_f[vals[6]], sLead_f[vals[5]], sLead_f[vals[4]], sLead_f[vals[3]], sLead_f[vals[2]], sLead_f[vals[1]], sLead_f[vals[0]]);
	St  = opCode(set_ps, sTrail_f[vals[7]], sTrail_f[vals[6]], sTrail_f[vals[5]], sTrail_f[vals[4]], sTrail_f[vals[3]], sTrail_f[vals[2]], sTrail_f[vals[1]], sTrail_f[vals[0]]);
	S   = opCode(add_ps, Sl, St);
	msk = opCode(cmp_ps, opCode(castsi256_ps, M), opCode(setzero_ps), _CMP_EQ_OQ);

	return	opCode(blendv_ps,
			opCode(mul_ps,
				opCode(cvtepi32_ps, opCode(sllv_epi32, two, opCode(sub_epi32, M, opCode(set1_epi32, 1)))),
				opCode(add_ps, Sl, opCode(add_ps, St, opCode(mul_ps, S, Q)))),
			opCode(add_ps, Sl, opCode(add_ps, St, opCode(mul_ps, S, Q))),
			msk);
}

#ifdef	__AVX512F__
inline _MData_	opCode(exp_ps, _MData_ &x)
{
	return	opCode(set_ps, std::exp(x[15]), std::exp(x[14]), std::exp(x[13]), std::exp(x[12]), std::exp(x[11]), std::exp(x[10]), std::exp(x[9]), std::exp(x[8]),
			       std::exp(x[7]),  std::exp(x[6]),  std::exp(x[5]),  std::exp(x[4]),  std::exp(x[3]),  std::exp(x[2]),  std::exp(x[1]), std::exp(x[0]));
}
#elif   defined(__AVX__)
inline _MData_	opCode(exp_ps, _MData_ &x)
{
	return	opCode(set_ps, std::exp(x[7]), std::exp(x[6]), std::exp(x[5]), std::exp(x[4]), std::exp(x[3]), std::exp(x[2]), std::exp(x[1]), std::exp(x[0]));
}
#else
inline _MData_	opCode(exp_ps, _MData_ &x)
{
	return	opCode(set_ps, std::exp(x[3]), std::exp(x[2]), std::exp(x[1]), std::exp(x[0]));
}
#endif

#endif

#undef	_MData_
#undef	_PREFIX
#undef opCode_P
#undef opCode_N
#undef opCode

#endif
