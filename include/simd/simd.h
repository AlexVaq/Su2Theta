#ifndef	__SIMD
	#define	__SIMD

	namespace Simd {

		#if   defined(AVX512F)
		static	size_t sAlign = 64;
		#elif defined(AVX)
		static	size_t sAlign = 32;
		#else
		static	size_t sAlign = 16;
		#endif

//		class	vSimd {};
	}
#endif
