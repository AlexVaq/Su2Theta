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

int	main (int argc, char *argv[]) {

	constexpr size_t nIters = 4294967296/128;
	constexpr size_t pFreq  = 268435456/64;

	initSu2 (argc, argv);

	initRandom();

	std::chrono::high_resolution_clock::time_point start, stop;
	std::chrono::nanoseconds avxElapsed, stdElapsed, othElapsed;

	start = std::chrono::high_resolution_clock::now();

	float x0 = genRand();
	Simd_f testVar4(x0, genRand(), genRand(), genRand(), genRand(), genRand(), genRand(), genRand());
	Simd_f zero    (0.);
	Simd_f testVar5(zero);
	Simd_f testVar (zero);

	_MData_ tVar;
	_MData_ tVar4;
	_MData_ tVar5;

	tVar  = opCode(set1_ps, 0.);
	tVar5 = opCode(set1_ps, 0.);
	tVar4 = opCode(set_ps,  testVar4[7],  testVar4[6],  testVar4[5],  testVar4[4],  testVar4[3],  testVar4[2],  testVar4[1],  testVar4[0]);

	for (size_t k=0; k<nIters; k++) {
		_MData_ tVar3 = tVar4;
		size_t i = k%16384;
		_MData_ tVar2 = opCode(set_ps, (float) i*0.98e-9f, (float) i*0.34e-8f, (float) i*0.49e-8f, (float) i*0.1e-8f, (float) i*0.14e-8f, (float) i *0.25e-8f, (float) i*0.75e-8f, (float) i*1.e-8f);

		tVar  = opCode(add_ps, tVar3, tVar5);
		tVar3 = opCode(mul_ps, tVar,  tVar2);
		tVar  = opCode(sub_ps, tVar3, tVar2);
		tVar  = opCode(sub_ps, tVar,  tVar5);
		tVar5 = tVar;
	}

	stop  = std::chrono::high_resolution_clock::now();

	othElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

	float dataOut[8];

	opCode(store_ps, dataOut, tVar5);

	printf ("Avx took %lu nanoseconds (intrinsics)\n", othElapsed);
	printf ("Result: %f %f %f %f %f %f %f %f\n", dataOut[0], dataOut[1], dataOut[2], dataOut[3], dataOut[4], dataOut[5], dataOut[6], dataOut[7]);

	start = std::chrono::high_resolution_clock::now();

	for (size_t k=0; k<nIters; k++) {
		Simd_f testVar3(testVar4);
		size_t i = k%16384;
		Simd_f testVar2((float) i*1e-8f, (float) i*0.75e-8f, (float) i*0.25e-8f, (float) i*0.14e-8f, (float) i*0.1e-8f, (float) i *0.49e-8f, (float) i*0.34e-8f, (float) i*0.98e-9f);

		for (int j=0; j<8; j++) {
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
	printf ("Result: %f %f %f %f %f %f %f %f\n", testVar5[0], testVar5[1], testVar5[2], testVar5[3], testVar5[4], testVar5[5], testVar5[6], testVar5[7]);

	testVar  = zero;
	testVar5 = zero;

	start = std::chrono::high_resolution_clock::now();

	for (size_t k=0; k<nIters; k++) {
		Simd_f testVar3(testVar4);
		size_t i = k%16384;
		Simd_f testVar2((float) i*1e-8f, (float) i*0.75e-8f, (float) i*0.25e-8f, (float) i*0.14e-8f, (float) i*0.1e-8f, (float) i *0.49e-8f, (float) i*0.34e-8f, (float) i*0.98e-9f);

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
}
