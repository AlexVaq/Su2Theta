#include <memory>
#include <chrono>
#include "enumFields.h"
#include "lattice/lattice.h"
#include "su2/su2.h"
#include "utils/memAlloc.h"
#include "utils/logger.h"
#include "utils/profiler.h"
#include "utils/misc.h"
#include "simd/simd.h"
#include "random/random.h"

using namespace Su2Rand;
using namespace Simd;

#ifdef	__AVX512F__
constexpr int nSimd = 16;
#else
constexpr int nSimd =  8;
#endif

typedef union ieee754_f {
	float  f;
	int    u;
}	ieee754_f;

typedef union ieee754_d {
	double f;
	size_t u;
}	ieee754_d;

int	main (int argc, char *argv[]) {

	constexpr size_t nIters = 4294967296/128;
	constexpr size_t pFreq  = 268435456/64;

	initSu2 (argc, argv);

	initRandom();

	std::chrono::high_resolution_clock::time_point start, stop;
	std::chrono::nanoseconds avxElapsed, stdElapsed, othElapsed;

	start = std::chrono::high_resolution_clock::now();

	float x0 = genRand();
#if	defined(__AVX512F__)
	Simd_f testVar4(x0, genRand(), genRand(), genRand(), genRand(), genRand(), genRand(), genRand(),
			x0, genRand(), genRand(), genRand(), genRand(), genRand(), genRand(), genRand());
#else
	Simd_f testVar4(x0, genRand(), genRand(), genRand(), genRand(), genRand(), genRand(), genRand());
#endif
	Simd_f zero    (0.);
	Simd_f testVar5(zero);
	Simd_f testVar (zero);

	_MData_ tVar;
	_MData_ tVar4;
	_MData_ tVar5;

	tVar  = opCode(set1_ps, 0.);
	tVar5 = opCode(set1_ps, 0.);
#if	defined(__AVX512F__)
	tVar4 = opCode(set_ps,  testVar4[15],  testVar4[14],  testVar4[13],  testVar4[12],  testVar4[11],  testVar4[10],  testVar4[9],  testVar4[8],
				testVar4[7],   testVar4[6],   testVar4[5],   testVar4[4],   testVar4[3],   testVar4[2],   testVar4[1],  testVar4[0]);
#else
	tVar4 = opCode(set_ps,  testVar4[7],  testVar4[6],  testVar4[5],  testVar4[4],  testVar4[3],  testVar4[2],  testVar4[1],  testVar4[0]);
#endif

	for (size_t k=0; k<nIters; k++) {
		_MData_ tVar3 = tVar4;
		size_t i = k%16384;
#if	defined(__AVX512F__)
		_MData_ tVar2 = opCode(set_ps, (float) i*0.98e-9f, (float) i*0.34e-8f, (float) i*0.49e-8f, (float) i*0.1e-8f, (float) i*0.14e-8f, (float) i *0.25e-8f, (float) i*0.75e-8f, (float) i*1.e-8f,
					       (float) i*0.98e-9f, (float) i*0.34e-8f, (float) i*0.49e-8f, (float) i*0.1e-8f, (float) i*0.14e-8f, (float) i *0.25e-8f, (float) i*0.75e-8f, (float) i*1.e-8f);
#else
		_MData_ tVar2 = opCode(set_ps, (float) i*0.98e-9f, (float) i*0.34e-8f, (float) i*0.49e-8f, (float) i*0.1e-8f, (float) i*0.14e-8f, (float) i *0.25e-8f, (float) i*0.75e-8f, (float) i*1.e-8f);
#endif

		tVar  = opCode(add_ps, tVar3, tVar5);
		tVar3 = opCode(mul_ps, tVar,  tVar2);
		tVar  = opCode(sub_ps, tVar3, tVar2);
		tVar  = opCode(sub_ps, tVar,  tVar5);
		tVar5 = tVar;
	}

	stop  = std::chrono::high_resolution_clock::now();

	othElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

	float dataOut[nSimd];

	opCode(store_ps, dataOut, tVar5);

	printf ("Avx took %lu nanoseconds (intrinsics)\n", othElapsed);
	printf ("Result: %f %f %f %f %f %f %f %f\n", dataOut[0], dataOut[1], dataOut[2],  dataOut[3],  dataOut[4],  dataOut[5],  dataOut[6],  dataOut[7]);
#if	defined(__AVX512F__)
	printf ("        %f %f %f %f %f %f %f %f\n", dataOut[8], dataOut[9], dataOut[10], dataOut[11], dataOut[12], dataOut[13], dataOut[14], dataOut[15]);
#endif

	start = std::chrono::high_resolution_clock::now();

	for (size_t k=0; k<nIters; k++) {
		Simd_f testVar3(testVar4);
		size_t i = k%16384;
#if	defined(__AVX512F__)
		Simd_f testVar2((float) i*1e-8f, (float) i*0.75e-8f, (float) i*0.25e-8f, (float) i*0.14e-8f, (float) i*0.1e-8f, (float) i *0.49e-8f, (float) i*0.34e-8f, (float) i*0.98e-9f,
				(float) i*1e-8f, (float) i*0.75e-8f, (float) i*0.25e-8f, (float) i*0.14e-8f, (float) i*0.1e-8f, (float) i *0.49e-8f, (float) i*0.34e-8f, (float) i*0.98e-9f);
#else
		Simd_f testVar2((float) i*1e-8f, (float) i*0.75e-8f, (float) i*0.25e-8f, (float) i*0.14e-8f, (float) i*0.1e-8f, (float) i *0.49e-8f, (float) i*0.34e-8f, (float) i*0.98e-9f);
#endif

		for (int j=0; j<nSimd; j++) {
			testVar[j]   = (float) testVar3[j]+ (float) testVar5[j];
			testVar3[j]  = (float) testVar[j] * (float) testVar2[j];
			testVar[j]   = (float) testVar3[j]- (float) testVar2[j];
			testVar[j]  -= (float) testVar5[j];
			testVar5[j]  = (float) testVar[j];
		}
	}

	stop  = std::chrono::high_resolution_clock::now();

	stdElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

	printf ("Standard took %lu nanoseconds\n", stdElapsed);
	printf ("Result: %f %f %f %f %f %f %f %f\n", testVar5[0], testVar5[1], testVar5[2],  testVar5[3],  testVar5[4],  testVar5[5],  testVar5[6],  testVar5[7]);
#if	defined(__AVX512F__)
	printf ("        %f %f %f %f %f %f %f %f\n", testVar5[8], testVar5[9], testVar5[10], testVar5[11], testVar5[12], testVar5[13], testVar5[14], testVar5[15]);
#endif

	testVar  = zero;
	testVar5 = zero;

	start = std::chrono::high_resolution_clock::now();

	for (size_t k=0; k<nIters; k++) {
		Simd_f testVar3(testVar4);
		size_t i = k%16384;
#if	defined(__AVX512F__)
		Simd_f testVar2((float) i*1e-8f, (float) i*0.75e-8f, (float) i*0.25e-8f, (float) i*0.14e-8f, (float) i*0.1e-8f, (float) i *0.49e-8f, (float) i*0.34e-8f, (float) i*0.98e-9f,
				(float) i*1e-8f, (float) i*0.75e-8f, (float) i*0.25e-8f, (float) i*0.14e-8f, (float) i*0.1e-8f, (float) i *0.49e-8f, (float) i*0.34e-8f, (float) i*0.98e-9f);
#else
		Simd_f testVar2((float) i*1e-8f, (float) i*0.75e-8f, (float) i*0.25e-8f, (float) i*0.14e-8f, (float) i*0.1e-8f, (float) i *0.49e-8f, (float) i*0.34e-8f, (float) i*0.98e-9f);
#endif
		testVar  = testVar3+ testVar5;
		testVar3 = testVar * testVar2;
		testVar  = testVar3- testVar2;
		testVar -= testVar5;
		testVar5 = testVar;
	}

	stop  = std::chrono::high_resolution_clock::now();

	avxElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

	printf ("Avx took %lu nanoseconds (class)\n", avxElapsed);
	printf ("Result: %f %f %f %f %f %f %f %f\n", testVar5[0], testVar5[1], testVar5[2], testVar5[3], testVar5[4], testVar5[5], testVar5[6], testVar5[7]);
#if	defined(__AVX512F__)
	printf ("        %f %f %f %f %f %f %f %f\n", testVar5[8], testVar5[9], testVar5[10], testVar5[11], testVar5[12], testVar5[13], testVar5[14], testVar5[15]);
#endif

	float x1 = 0., x5 = 0.;

	for (size_t k=0; k<nIters; k++) {
		float x3 = x0;
		size_t i = k%16384;
		float x2 = (float) i*1e-8f;

		x1  = x3+x5;
		x3  = x1*x2;
		x1  = x3-x2;
		x1 -= x5;
		x5  = x1;
	}
	printf ("Result: %f\n", x5);
	printf ("Speedup\t\tx%.2lf\n", ((double) stdElapsed.count())/((double) avxElapsed.count()));
	printf ("Class vs intrinsics\tx%.2lf\n", ((double) othElapsed.count())/((double) avxElapsed.count()));

	printf("Logarithm test\n");

	std::vector<_MData_> base(nIters);
	std::vector<_MData_> news(nIters);
	std::vector<_MData_> outs(nIters);

	start = std::chrono::high_resolution_clock::now();

	for (size_t i=1; i<nIters; i++) {
#if	defined(__AVX512F__)
		base[i] = opCode(set_ps, std::log(((float) i)*7e-4f),   std::log(((float) i)*7e-5f),     std::log(((float) i)*0.25e-5f),  std::log(((float) i)*0.14e-6f),
					 std::log(((float) i)*7e-4f),   std::log(((float) i)*7e-5f),     std::log(((float) i)*0.25e-5f),  std::log(((float) i)*0.14e-6f),
					 std::log(((float) i)*0.1e-8f), std::log(((float) i)*0.49e-18f), std::log(((float) i)*0.34e-28f), std::log(((float) i)*0.98e-39f),
					 std::log(((float) i)*0.1e-8f), std::log(((float) i)*0.49e-18f), std::log(((float) i)*0.34e-28f), std::log(((float) i)*0.98e-39f));
#else
//		base[i] = opCode(set_ps, std::log(((float) i)*7e-4f),   std::log(((float) i)*7e-5f),     std::log(((float) i)*0.25e-5f),  std::log(((float) i)*0.14e-6f),
//					 std::log(((float) i)*0.1e-8f), std::log(((float) i)*0.49e-18f), std::log(((float) i)*0.34e-28f), std::log(((float) i)*0.98e-39f));
		base[i] = opCode(set_ps, std::log(((float) i)*7e-4f), std::log(((float) i)*7e-5f), std::log(((float) i)*0.34e-28f), std::log(((float) i)*0.98e-39f),
					 std::log(((float) i)*7e-4f), std::log(((float) i)*7e-5f), std::log(((float) i)*0.34e-28f), std::log(((float) i)*0.98e-39f));
#endif
	}

	stop  = std::chrono::high_resolution_clock::now();
	stdElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	printf ("Std took %lu nanoseconds\n", stdElapsed); fflush(stdout);

	start = std::chrono::high_resolution_clock::now();

	for (size_t i=1; i<nIters; i++) {
#if	defined(__AVX512F__)
		_MData_ tVar = opCode(set_ps, ((float) i)*7e-4f,   ((float) i)*7e-5f,     ((float) i)*0.25e-5f,  ((float) i)*0.14e-6f,
					      ((float) i)*7e-4f,   ((float) i)*7e-5f,     ((float) i)*0.25e-5f,  ((float) i)*0.14e-6f,
					      ((float) i)*0.1e-8f, ((float) i)*0.49e-18f, ((float) i)*0.34e-28f, ((float) i)*0.98e-39f,
					      ((float) i)*0.1e-8f, ((float) i)*0.49e-18f, ((float) i)*0.34e-28f, ((float) i)*0.98e-39f);
#else
		_MData_ tVar = opCode(set_ps, ((float) i)*7e-4f, ((float) i)*7e-5f, ((float) i)*0.34e-28f, ((float) i)*0.98e-39f,
					      ((float) i)*7e-4f, ((float) i)*7e-5f, ((float) i)*0.34e-28f, ((float) i)*0.98e-39f);
//		_MData_ tVar = opCode(set_ps, ((float) i)*7e-4f,   ((float) i)*7e-5f,     ((float) i)*0.25e-5f,  ((float) i)*0.14e-6f,
//					      ((float) i)*0.1e-8f, ((float) i)*0.49e-18f, ((float) i)*0.34e-28f, ((float) i)*0.98e-39f);
#endif
		news[i] = opCode(log_ps, tVar);
	}

	stop  = std::chrono::high_resolution_clock::now();
	avxElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

	printf ("Speedup\t\tx%.2lf\n", ((double) stdElapsed.count())/((double) avxElapsed.count()));
	printf ("Avx took %lu nanoseconds\n", avxElapsed); fflush(stdout);
	float	fMax = 0.f;
	int	iMax = 0;
	int	kMax = 0;
	for (size_t i=1; i<nIters; i++) {
	/*	if ((i%1048576) == 0) {
			printsVar(base[i], "base");
			printsVar(news[i], "news");
		}*/
		outs[i] = opCode(div_ps, opCode(sub_ps, base[i], news[i]), base[i]);
		for (int k=0; k<nSimd; k++)
			if (std::abs(outs[i][k]) > fMax) {
				fMax = std::abs(outs[i][k]);
				iMax = i;
				kMax = k;
			}
	}
	printf ("Max error (%d, %d) %.8e (M %.8e S %.8e)\n", iMax, kMax, fMax, news[iMax][kMax], base[iMax][kMax]); fflush(stdout);

#if	defined(__AVX512F__)
	std::vector<__m512d> based(nIters);
	std::vector<__m512d> newsd(nIters);
	std::vector<__m512d> outsd(nIters);
#else
	std::vector<__m256d> based(nIters);
	std::vector<__m256d> newsd(nIters);
	std::vector<__m256d> outsd(nIters);
#endif

	start = std::chrono::high_resolution_clock::now();

	for (size_t i=1; i<nIters; i++) {
#if	defined(__AVX512F__)
		based[i] = opCode(set_pd, std::log(((double) i)*7e-4),   std::log(((double) i)*7e-5),     std::log(((double) i)*0.25e-5),  std::log(((double) i)*0.14e-6),
					  std::log(((double) i)*0.1e-8), std::log(((double) i)*0.49e-18), std::log(((double) i)*0.34e-28), std::log(((double) i)*0.98e-39));
#else
		based[i] = opCode(set_pd, std::log(((double) i)*7e-4),   std::log(((double) i)*7e-5),     std::log(((double) i)*0.34e-28), std::log(((double) i)*0.98e-39));
#endif
	}

	stop  = std::chrono::high_resolution_clock::now();
	stdElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	printf ("Std took %lu nanoseconds\n", stdElapsed); fflush(stdout);

	start = std::chrono::high_resolution_clock::now();

	for (size_t i=1; i<nIters; i++) {
#if	defined(__AVX512F__)
		__m512d dVar = opCode(set_pd, ((double) i)*7e-4,   ((double) i)*7e-5,     ((double) i)*0.25e-5,  ((double) i)*0.14e-6,
					      ((double) i)*0.1e-8, ((double) i)*0.49e-18, ((double) i)*0.34e-28, ((double) i)*0.98e-39);
#else
		__m256d dVar = opCode(set_pd, ((double) i)*7e-4,   ((double) i)*7e-5,     ((double) i)*0.34e-28, ((double) i)*0.98e-39);
#endif
		newsd[i] = opCode(log_pd, dVar);
	}

	stop  = std::chrono::high_resolution_clock::now();
	avxElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

	printf ("Speedup\t\tx%.2lf\n", ((double) stdElapsed.count())/((double) avxElapsed.count()));
	printf ("Avx took %lu nanoseconds\n", avxElapsed); fflush(stdout);
	double	dMax = 0.0;
	for (size_t i=1; i<nIters; i++) {
	/*	if ((i%1048576) == 0) {
			printdVar(based[i], "base");
			printdVar(newsd[i], "news");
		}*/
		outsd[i] = opCode(div_pd, opCode(sub_pd, based[i], newsd[i]), based[i]);
		for (int k=0; k<nSimd/2; k++)
			if (std::abs(outs[i][k]) > dMax) {
				dMax = std::abs(outsd[i][k]);
				iMax = i;
				kMax = k;
			}
	}
	printf ("Max error %d %.16le (M %.16le S %.16le)\n", iMax, dMax, newsd[iMax][kMax], based[iMax][kMax]); fflush(stdout);
}
