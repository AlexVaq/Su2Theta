#include <memory>
#include "enumFields.h"
#include "lattice/lattice.h"
#include "su2/su2.h"
#include "su2/vSu2.h"
#include "action/action.h"
#include "simd/simd.h"
#include "utils/memAlloc.h"
#include "utils/logger.h"
#include "utils/profiler.h"
#include "utils/misc.h"
#include <chrono>

using namespace Simd;

constexpr int nTerm  = 16;//384;
constexpr int nIters = 32;//512;

int	main (int argc, char *argv[]) {
	std::chrono::high_resolution_clock::time_point start, stop;
	std::chrono::nanoseconds avxElapsed, stdElapsed, othElapsed;

	initSu2 (argc, argv);

	Lattice<Su2<float>,EvenOdd> *myLat = new Lattice<Su2<float>,EvenOdd>(16, 16);
	myLat->SetRand();

	Su2Action::Action   <Su2<float>,EvenOdd> wAction(*myLat, 2.05, 0.0);

	auto	sAct = wAction.allPts();

	printf("Plaquette value: %le\n\n", sAct/(6.*myLat->Volume()));

	auto	sMetro = Su2Action::Metropolis<Su2<float>,EvenOdd>(wAction);
	auto	wPlq   = Su2Action::Plaquette <Su2<float>,EvenOdd>(*myLat);

	Su2Tune::Tune(sMetro);
	Su2Tune::Tune(wPlq);

	start = std::chrono::high_resolution_clock::now();

	double	stdAvg = 0., stdErr = 0.;

	for	(int i = 0; i<nTerm; i++)
		sMetro();

	for	(int i = 0; i<nIters; i++) {
		sMetro();

		sAct = wPlq();

		stdAvg += sAct;
		stdErr += sAct*sAct;

		printf("%04d\t%.4f\n", i, sAct);
	}

	stop  = std::chrono::high_resolution_clock::now();

	delete myLat;

	stdAvg /= ((double) nIters);
	stdErr /= ((double) nIters);
	stdElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	printf("Elapsed time %.3lf ms\n", ((double) stdElapsed.count())/1e6);
	printf("Final Plaquette: %le +/- %le\n", stdAvg, std::sqrt((stdErr - stdAvg*stdAvg)/((double) (nIters - 1))));

	Lattice<vSu2<Simd_f>,EvenOdd> *myVLat = new Lattice<vSu2<Simd_f>,EvenOdd>(16, 16);
	myVLat->SetRand();

	Su2Action::Action<vSu2<Simd_f>,EvenOdd> wVAction(*myVLat, 2.05, 0.0);

	auto sVAct = wVAction.allPts();

	printf("Plaquette value: %le\n\n", sVAct/(6.*myVLat->Volume()));

	auto vMetro = Su2Action::Metropolis<vSu2<Simd_f>,EvenOdd>(wVAction);
	auto vPlq   = Su2Action::Plaquette <vSu2<Simd_f>,EvenOdd>(*myVLat);

	Su2Tune::Tune(vMetro);
	Su2Tune::Tune(vPlq);

	start = std::chrono::high_resolution_clock::now();

	double	avxAvg = 0., avxErr = 0.;

	for	(int i = 0; i<nTerm; i++)
		vMetro();

	for	(int i = 0; i<nIters; i++) {
		vMetro();

		sVAct = vPlq();

		avxAvg += sVAct;
		avxErr += sVAct*sVAct;

//		if (!(i%16)) {
      			printf("%04d\t%.4f\n", i, sVAct);
//		}

	}

	stop  = std::chrono::high_resolution_clock::now();

	delete myVLat;

	avxAvg /= ((double) nIters);
	avxErr /= ((double) nIters);
	avxElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	printf("Elapsed time %.3lf ms\n", ((double) avxElapsed.count())/1e6);
	printf("Final Plaquette: %le +/- %le\n", avxAvg, std::sqrt((avxErr - avxAvg*avxAvg)/((double) (nIters - 1))));
	printf("Difference %le\n", std::abs(avxAvg - stdAvg));
	printf("Speedup x%.2lf\n", ((double) stdElapsed.count())/((double) avxElapsed.count()));

	endSu2  ();
}
