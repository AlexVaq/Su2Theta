#include <simd/simd-tri.h>
#include <simd/simd-math.h>

#ifdef  __AVX512F__
	#include <simd/simd-Avx512.h>
#elif   defined(__AVX2__)
	#include <simd/simd-Avx2.h>
#elif   defined(__AVX__)
	#include <simd/simd-Avx.h>
#else
	#include <simd/simd-Sse4.h>
#endif

