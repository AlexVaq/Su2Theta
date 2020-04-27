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

#define	Float    double  // float
#define	vFloat   Simd_d  // Simd_f
#define Ordering Colored // EvenOdd

using namespace Simd;

constexpr int nTerm  = 0;
constexpr int nIters = 16;//512;
constexpr int nOvHB  = 3;

int	main (int argc, char *argv[]) {
	std::chrono::high_resolution_clock::time_point start, stop;
	std::chrono::nanoseconds avxElapsed, stdElapsed, othElapsed;

	initSu2 (argc, argv);

	auto *myLat = new Lattice<Su2<Float>,Ordering>(16, 16);
	myLat->SetRand();

	Su2Action::Action<Su2<Float>,Ordering> wAction(*myLat, 2.05, 0.0);

	printf("Start tuning\n"); fflush(stdout);
	Su2Tune::Tune(wAction);

	auto	sAct = wAction.allPts();
	printf("Plaquette value: %le\n\n", sAct/(6.*myLat->Volume()));

	double	stdAvg = 0., stdErr = 0.;

	auto	sHB = Su2Action::HeatBath <Su2<Float>,Ordering>(wAction);
	auto	sOR = Su2Action::OverRelax<Su2<Float>,Ordering>(wAction);
	auto	sPq = Su2Action::Plaquette<Su2<Float>,Ordering>(*myLat);

	Su2Tune::Tune(sHB);
	Su2Tune::Tune(sOR);
	Su2Tune::Tune(sPq);

	sAct = sPq();
	printf("Plaquette before: %le\n\n", sAct);
	start = std::chrono::high_resolution_clock::now();

	for	(int i = 0; i<nTerm; i++) {
		sHB(nOvHB);
		sOR();
	}

	for	(int i = 0; i<nIters; i++) {
//		sHB(nOvHB);
		sOR();

		sAct = sPq();

		stdAvg += sAct;
		stdErr += sAct*sAct;

//		if (!(i%16)) {
			printf("%04d\t%.7e\n", i, sAct);
//		}

	}

	stop  = std::chrono::high_resolution_clock::now();
	sAct = sPq();
	printf("Plaquette after:  %le\n\n", sAct);

	delete myLat;

	stdAvg /= ((double) nIters);
	stdErr /= ((double) nIters);
	stdElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	printf("Elapsed time %.3lf ms\n", ((double) stdElapsed.count())/1e6);
	printf("Final Plaquette: %le +/- %le\n", stdAvg, std::sqrt((stdErr - stdAvg*stdAvg)/((double) (nIters - 1))));

	Lattice<vSu2<vFloat>,Ordering> *myVLat = new Lattice<vSu2<vFloat>,Ordering>(16, 16);
	myVLat->SetRand();

	Su2Action::Action<vSu2<vFloat>,Ordering> wVAction(*myVLat, 2.05, 0.0);

	printf("Start tuning\n"); fflush(stdout);
	Su2Tune::Tune(wVAction);
	auto sVAct = wVAction.allPts();

	printf("Plaquette value: %le\n\n", sVAct/(6.*myVLat->Volume()));

	double	avxAvg = 0., avxErr = 0.;

	auto	vHB = Su2Action::HeatBath <vSu2<vFloat>,Ordering>(wVAction);
	auto	vOR = Su2Action::OverRelax<vSu2<vFloat>,Ordering>(wVAction);
	auto	vPq = Su2Action::Plaquette<vSu2<vFloat>,Ordering>(*myVLat);

	Su2Tune::Tune(vHB);
	Su2Tune::Tune(vOR);
	Su2Tune::Tune(vPq);

	start = std::chrono::high_resolution_clock::now();

	for	(int i = 0; i<nTerm; i++) {
		vHB(nOvHB);
		vOR();
	}

	for	(int i = 0; i<nIters; i++) {
//		vHB(nOvHB);
		vOR();

		sVAct = vPq();

		avxAvg += sVAct;
		avxErr += sVAct*sVAct;

	//	if (!(i%16)) {
			printf("%04d\t%.7e\n", i, sVAct);
	//	}

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
