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

	constexpr size_t nIters = 4294967296/64;
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

	printf ("Std took %lu nanoseconds\n", stdElapsed);
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
	printf ("Class vs intr\tx%.2lf\n", ((double) othElapsed.count())/((double) avxElapsed.count()));

	printf("\n\nLogarithm test\n\n");
	printf("Single precision\n");

	double rangeLowf  = 1e-43;	// Min 1.4e-45;
	double rangeHighf = 3e+38;	// Max 3.4e+38;
	double stepf	  = std::exp((std::log(rangeHighf) - std::log(rangeLowf))/((double) nIters));
	double cValf	  = rangeLowf;

	void  *insf;
	void  *news;
	std::vector<float> base(nIters);
	posix_memalign(&insf, Simd::sAlign, nIters*sizeof(float));
	posix_memalign(&news, Simd::sAlign, nIters*sizeof(float));

	for (int i=0; i<nIters; i++) {
		static_cast<float*>(insf)[i] = cValf;
		cValf *= stepf;
	}

	start = std::chrono::high_resolution_clock::now();

	for (size_t i=0; i<nIters; i++)
		base[i] = std::log(static_cast<float*>(insf)[i]);

	stop  = std::chrono::high_resolution_clock::now();
	stdElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	printf ("Std took %lu nanoseconds\n", stdElapsed); fflush(stdout);

	start = std::chrono::high_resolution_clock::now();

	for (size_t i=0; i<nIters; i+=Simd_f::sWide)
		opCode(store_ps, &(static_cast<float*>(news)[i]), opCode(log_ps, opCode(load_ps, &(static_cast<float*>(insf)[i]))));

	stop  = std::chrono::high_resolution_clock::now();
	avxElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

	printf ("Avx took %lu nanoseconds\n", avxElapsed); fflush(stdout);
	printf ("Speedup\t\tx%.2lf\n", ((double) stdElapsed.count())/((double) avxElapsed.count()));

	double	sMaxf = 0.f;
	double	aMaxf = 0.f;
	size_t	iaMax = 0;
	size_t	isMax = 0;

	for (size_t i=0; i<nIters; i++) {
		double tLog = std::log((double) static_cast<float*>(insf)[i]);
		double sTmp = std::abs((((double) static_cast<float*>(news)[i]) - tLog)/tLog);
		double aTmp = std::abs((((double) base[i]) - tLog)/tLog);
		if (sTmp > sMaxf) {
			sMaxf = sTmp;
			isMax = i;
		}
		if (aTmp > aMaxf) {
			aMaxf = aTmp;
			iaMax = i;
		}
	}
	printf ("Max Avx error %.9e @ %.9e C %.12le T %.12le\n", sMaxf, static_cast<float*>(insf)[isMax], static_cast<float*>(news)[isMax],
								 std::log((double) static_cast<float*>(insf)[isMax]));
	printf ("Max Std error %.9e @ %.9e C %.12le T %.12le\n", aMaxf, static_cast<float*>(insf)[iaMax], base[iaMax], std::log((double) static_cast<float*>(insf)[iaMax]));
	printf ("Range [%.9e %.9e - %.9e %.9e] - %.9e\n\n", static_cast<float*>(insf)[0],        static_cast<float*>(insf)[1],
							    static_cast<float*>(insf)[nIters-2], static_cast<float*>(insf)[nIters-1], stepf);
       	fflush(stdout);

	printf("\nDouble precision\n");

	void	*insd;
	void	*newsd;
	std::vector<double> based(nIters);
	std::vector<double> outsd(nIters);
	posix_memalign(&insd,  Simd::sAlign, nIters*sizeof(double));
	posix_memalign(&newsd, Simd::sAlign, nIters*sizeof(double));

	long double rangeLowd  = 2e-307;	// Min 4.9407e-324 but just in case
	long double rangeHighd = 1e+308;	// Max 3.59e+308 but it won't take it --> Inf
	long double stepd      = std::exp((std::log((long double) rangeHighd) - std::log((long double) rangeLowd))/((long double) nIters));
	long double cVald      = rangeLowd;

	for (int i=0; i<nIters; i++) {
		static_cast<double*>(insd)[i] = (double) cVald;
		cVald *= stepd;
	}

	start = std::chrono::high_resolution_clock::now();

	for (size_t i=0; i<nIters; i++)
		base[i] = std::log(static_cast<double*>(insd)[i]);

	stop  = std::chrono::high_resolution_clock::now();
	stdElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	printf ("Std took %lu nanoseconds\n", stdElapsed); fflush(stdout);

	start = std::chrono::high_resolution_clock::now();

	for (size_t i=0; i<nIters; i+=Simd_d::sWide)
		opCode(store_pd, &(static_cast<double*>(newsd)[i]), opCode(log_pd, opCode(load_pd, &(static_cast<double*>(insd)[i]))));

	stop  = std::chrono::high_resolution_clock::now();
	avxElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

	printf ("Avx took %lu nanoseconds\n", avxElapsed); fflush(stdout);
	printf ("Speedup\t\tx%.2lf\n", ((double) stdElapsed.count())/((double) avxElapsed.count()));
	double	sMaxd = 0.0;
	double	aMaxd = 0.0;
	for (size_t i=0; i<nIters; i++) {
		long double tLog = std::log(  (long double) (static_cast<double*>(insd)[i]));
		long double sTmp = std::abs((((long double) (static_cast<double*>(newsd)[i])) - tLog)/tLog);
		long double aTmp = std::abs((((long double) (based[i])) - tLog)/tLog);
		if (sTmp > sMaxd) {
			sMaxd = sTmp;
			isMax = i;
		}
		if (aTmp > aMaxd) {
			aMaxd = aTmp;
			iaMax = i;
		}
	}
	printf ("Max Avx error %.12le @ %.12le C %.12le T %.16Le\n", sMaxd, static_cast<double*>(insd)[isMax], static_cast<double*>(newsd)[isMax],
								     std::log((long double) (static_cast<double*>(insd)[isMax])));
	printf ("Max Std error %.12le @ %.12le C %.12le T %.16Le\n", aMaxd, static_cast<double*>(insd)[iaMax], based[iaMax], std::log((long double) (static_cast<double*>(insd)[iaMax])));
	printf ("Range [%.12le %.12le - %.12le %.12le] - %.16Le\n\n", static_cast<double*>(insd)[0],        static_cast<double*>(insd)[1],
								      static_cast<double*>(insd)[nIters-2], static_cast<double*>(insd)[nIters-1], stepd);
       	fflush(stdout);

	printf("\n\nExponential test\n\n");
	printf("Single precision\n");

	for (size_t i=0; i<nIters; i++)
		base[i] = std::exp(static_cast<float*>(insf)[i]);

	stop  = std::chrono::high_resolution_clock::now();
	stdElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	printf ("Std took %lu nanoseconds\n", stdElapsed); fflush(stdout);

	start = std::chrono::high_resolution_clock::now();

	for (size_t i=0; i<nIters; i+=Simd_f::sWide)
		opCode(store_ps, &(static_cast<float*>(news)[i]), opCode(exp_ps, opCode(load_ps, &(static_cast<float*>(insf)[i]))));

	stop  = std::chrono::high_resolution_clock::now();
	avxElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

	printf ("Avx took %lu nanoseconds\n", avxElapsed); fflush(stdout);
	printf ("Speedup\t\tx%.2lf\n", ((double) stdElapsed.count())/((double) avxElapsed.count()));

	sMaxf = 0.f;
	aMaxf = 0.f;

	for (size_t i=0; i<nIters; i++) {
		double tExp = std::exp((double) static_cast<float*>(insf)[i]);

		if (tExp > rangeHighf)
			continue;

		double sTmp = std::abs((((double) static_cast<float*>(news)[i]) - tExp)/tExp);
		double aTmp = std::abs((((double) base[i]) - tExp)/tExp);
		if (sTmp > sMaxf) {
			sMaxf = sTmp;
			isMax = i;
		}
		if (aTmp > aMaxf) {
			aMaxf = aTmp;
			iaMax = i;
		}
	}
	printf ("Max Avx error %.9e @ %.9e C %.12le T %.12le\n", sMaxf, static_cast<float*>(insf)[isMax], static_cast<float*>(news)[isMax],
								 std::exp((double) static_cast<float*>(insf)[isMax]));
	printf ("Max Std error %.9e @ %.9e C %.12le T %.12le\n", aMaxf, static_cast<float*>(insf)[iaMax], base[iaMax], std::exp((double) static_cast<float*>(insf)[iaMax]));
	printf ("Range [%.9e %.9e - %.9e %.9e] - %.9e\n\n", static_cast<float*>(insf)[0],        static_cast<float*>(insf)[1],
							    static_cast<float*>(insf)[nIters-2], static_cast<float*>(insf)[nIters-1], stepf);
       	fflush(stdout);

	for (size_t i=0; i<nIters; i++)
		static_cast<float*>(insf)[i] *= (-1.0);

	for (size_t i=0; i<nIters; i++)
		base[i] = std::exp(static_cast<float*>(insf)[i]);

	stop  = std::chrono::high_resolution_clock::now();
	stdElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	printf ("Std took %lu nanoseconds\n", stdElapsed); fflush(stdout);

	start = std::chrono::high_resolution_clock::now();

	for (size_t i=0; i<nIters; i+=Simd_f::sWide)
		opCode(store_ps, &(static_cast<float*>(news)[i]), opCode(exp_ps, opCode(load_ps, &(static_cast<float*>(insf)[i]))));

	stop  = std::chrono::high_resolution_clock::now();
	avxElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

	printf ("Avx took %lu nanoseconds\n", avxElapsed); fflush(stdout);
	printf ("Speedup\t\tx%.2lf\n", ((double) stdElapsed.count())/((double) avxElapsed.count()));

	sMaxf = 0.f;
	aMaxf = 0.f;

	for (size_t i=0; i<nIters; i++) {
		double tExp = std::exp((double) static_cast<float*>(insf)[i]);
		double sTmp = std::abs((((double) static_cast<float*>(news)[i]) - tExp)/tExp);
		double aTmp = std::abs((((double) base[i]) - tExp)/tExp);

		if (tExp < rangeLowf)
			continue;

		if (sTmp > sMaxf) {
			sMaxf = sTmp;
			isMax = i;
		}
		if (aTmp > aMaxf) {
			aMaxf = aTmp;
			iaMax = i;
		}
	}
	printf ("Max Avx error %.9e @ %.9e C %.12le T %.12le\n", sMaxf, static_cast<float*>(insf)[isMax], static_cast<float*>(news)[isMax],
								 std::exp((double) static_cast<float*>(insf)[isMax]));
	printf ("Max Std error %.9e @ %.9e C %.12le T %.12le\n", aMaxf, static_cast<float*>(insf)[iaMax], base[iaMax], std::exp((double) static_cast<float*>(insf)[iaMax]));
	printf ("Range [%.9e %.9e - %.9e %.9e] - %.9e\n\n", static_cast<float*>(insf)[0],        static_cast<float*>(insf)[1],
							    static_cast<float*>(insf)[nIters-2], static_cast<float*>(insf)[nIters-1], stepf);
       	fflush(stdout);

	printf("\nDouble precision\n");

	for (size_t i=0; i<nIters; i++)
		based[i] = std::exp(static_cast<double*>(insd)[i]);

	stop  = std::chrono::high_resolution_clock::now();
	stdElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	printf ("Std took %lu nanoseconds\n", stdElapsed); fflush(stdout);

	start = std::chrono::high_resolution_clock::now();

	for (size_t i=0; i<nIters; i+=Simd_d::sWide)
		opCode(store_pd, &(static_cast<double*>(newsd)[i]), opCode(exp_pd, opCode(load_pd, &(static_cast<double*>(insd)[i]))));

	stop  = std::chrono::high_resolution_clock::now();
	avxElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

	printf ("Avx took %lu nanoseconds\n", avxElapsed); fflush(stdout);
	printf ("Speedup\t\tx%.2lf\n", ((double) stdElapsed.count())/((double) avxElapsed.count()));

	sMaxd = 0.0;
	aMaxd = 0.0;

	for (size_t i=0; i<nIters; i++) {
		long double tExp = std::exp((long double) (static_cast<double*>(insd)[i]));
		long double sTmp = std::abs((((long double) (static_cast<double*>(newsd)[i])) - tExp)/tExp);
		long double aTmp = std::abs((((long double) (based[i])) - tExp)/tExp);

		if (tExp > rangeHighd)
			continue;

		if (sTmp > sMaxd) {
			sMaxd = sTmp;
			isMax = i;
		}
		if (aTmp > aMaxd) {
			aMaxd = aTmp;
			iaMax = i;
		}
	}
	printf ("Max Avx error %.12le @ %.12le C %.12le T %.16Le\n", sMaxd, static_cast<double*>(insd)[isMax], static_cast<double*>(newsd)[isMax],
								     std::exp((long double) (static_cast<double*>(insd)[isMax])));
	printf ("Max Std error %.12le @ %.12le C %.12le T %.16Le\n", aMaxd, static_cast<double*>(insd)[iaMax], based[iaMax], std::exp((long double) (static_cast<double*>(insd)[iaMax])));
	printf ("Range [%.12le %.12le - %.12le %.12le] - %.16Le\n\n", static_cast<double*>(insd)[0],        static_cast<double*>(insd)[1],
								      static_cast<double*>(insd)[nIters-2], static_cast<double*>(insd)[nIters-1], stepd);
       	fflush(stdout);

	for (size_t i=0; i<nIters; i++)
		static_cast<double*>(insd)[i] *= (-1.0);

	for (size_t i=0; i<nIters; i++)
		based[i] = std::exp(static_cast<double*>(insd)[i]);

	stop  = std::chrono::high_resolution_clock::now();
	stdElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	printf ("Std took %lu nanoseconds\n", stdElapsed); fflush(stdout);

	start = std::chrono::high_resolution_clock::now();

	for (size_t i=0; i<nIters; i+=Simd_d::sWide)
		opCode(store_pd, &(static_cast<double*>(newsd)[i]), opCode(exp_pd, opCode(load_pd, &(static_cast<double*>(insd)[i]))));

	stop  = std::chrono::high_resolution_clock::now();
	avxElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

	printf ("Avx took %lu nanoseconds\n", avxElapsed); fflush(stdout);
	printf ("Speedup\t\tx%.2lf\n", ((double) stdElapsed.count())/((double) avxElapsed.count()));

	sMaxd = 0.0;
	aMaxd = 0.0;

	for (size_t i=0; i<nIters; i++) {
		long double tExp = std::exp(  (long double) (static_cast<double*>(insd)[i]));
		long double sTmp = std::abs((((long double) (static_cast<double*>(newsd)[i])) - tExp)/tExp);
		long double aTmp = std::abs((((long double) (based[i])) - tExp)/tExp);

		if (tExp < rangeLowd)
			continue;

		if (sTmp > sMaxd) {
			sMaxd = sTmp;
			isMax = i;
		}
		if (aTmp > aMaxd) {
			aMaxd = aTmp;
			iaMax = i;
		}
	}

	printf ("Max Avx error %.12le @ %.12le C %.12le T %.16Le\n", sMaxd, static_cast<double*>(insd)[isMax], static_cast<double*>(newsd)[isMax],
								     std::exp((long double) (static_cast<double*>(insd)[isMax])));
	printf ("Max Std error %.12le @ %.12le C %.12le T %.16Le\n", aMaxd, static_cast<double*>(insd)[iaMax], based[iaMax], std::exp((long double) (static_cast<double*>(insd)[iaMax])));
	printf ("Range [%.12le %.12le - %.12le %.12le] - %.16Le\n\n", static_cast<double*>(insd)[0],        static_cast<double*>(insd)[1],
								      static_cast<double*>(insd)[nIters-2], static_cast<double*>(insd)[nIters-1], stepd);
       	fflush(stdout);

	printf ("\nRandom number generator test\n"); fflush(stdout);
	for (size_t i=0; i<16; i++) {
		printsVar((genVRand<Simd_f>()).raw(), "Random");
	}

	(exp(Simd_d(-700.0))).Print("-700");
	(exp(Simd_d(-708.0))).Print("-708");
	(exp(Simd_d(-708.4))).Print("-708.4");
	(exp(Simd_d(-708.5))).Print("-708.5");
	(exp(Simd_d(-710.0))).Print("-710.0");
	(exp(Simd_d(-715.0))).Print("-715.0");

	free(news); free(newsd);
	free(insf); free(insd);

	for (int i=0; i<nIters; i++) {
		auto rnd = vGenRand<Simd_d>();
		rnd.Print("Rng: ");
	}
}
