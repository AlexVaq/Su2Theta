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

constexpr int nTerm  = 4;
constexpr int nIters = 256;

int	main (int argc, char *argv[]) {
	std::chrono::high_resolution_clock::time_point start, stop;
	std::chrono::nanoseconds avxElapsed, stdElapsed, othElapsed;

	initSu2 (argc, argv);

	Lattice<Su2<double>> *myLat = new Lattice<Su2<double>>(32, 32, PrecSingle);
	myLat->SetRand();

	Su2Action::Action<Su2<double>> wAction(*myLat, 6.0, 0.0, false);

	auto	sAct = wAction.allPts();

	printf("Plaquette value: %le\nPoints %zu (%zu)\n", sAct/myLat->Volume(), myLat->Volume(), myLat->oVol());

	start = std::chrono::high_resolution_clock::now();

	double	stdAvg = 0., stdErr = 0.;

	for	(int i = 0; i<nTerm; i++) {
		Su2HB::HeatBath(wAction);
	}

	for	(int i = 0; i<nIters; i++) {
		Su2HB::HeatBath(wAction);
		sAct = wAction.allPts()/myLat->Volume();

		stdAvg += sAct;
		stdErr += sAct*sAct;

		if (!(i%16)) {
			printf("%04d\t%.4f\n", i, sAct);
		}
	}

	stop  = std::chrono::high_resolution_clock::now();

	delete myLat;

	stdAvg /= ((double) nIters);
	stdErr /= ((double) nIters);
	stdElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	printf("Elapsed time %.3lf ms\n", ((double) stdElapsed.count())/1e6);
	printf("Final Plaquette: %le +/- %le\n", stdAvg, std::sqrt((stdErr - stdAvg*stdAvg)/((double) (nIters - 1))));

	Lattice<vSu2<Simd_d>> *myVLat = new Lattice<vSu2<Simd_d>>(32, 32, PrecSingle);
	myVLat->SetRand();

	Su2Action::Action<vSu2<Simd_d>> wVAction(*myVLat, 6.0, 0.0, false);

	auto sVAct = wVAction.allPts();

	printf("Plaquette value: %le\nPoints %zu (%zu)\n", sVAct/myVLat->Volume(), myVLat->Volume(), myVLat->oVol());

	start = std::chrono::high_resolution_clock::now();

	double	avxAvg = 0., avxErr = 0.;

	for	(int i = 0; i<nTerm; i++) {
		Su2HB::HeatBath(wVAction);
	}

	for	(int i = 0; i<nIters; i++) {
		Su2HB::HeatBath(wVAction);
		sVAct = wVAction.allPts()/myVLat->Volume();

		avxAvg += sVAct;
		avxErr += sVAct*sVAct;

		if (!(i%16)) {
			printf("%04d\t%.4f\n", i, sVAct);
		}
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
