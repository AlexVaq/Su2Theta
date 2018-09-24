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

constexpr int nTerm  = 1024;
constexpr int nIters = 32768;
constexpr int nOvHB  = 3;
constexpr int nTries = 4;

int	main (int argc, char *argv[]) {
	std::chrono::high_resolution_clock::time_point start, stop;
	std::chrono::nanoseconds avxElapsed, stdElapsed, othElapsed;

	initSu2 (argc, argv);

	LogOut("Test plaquette action\n");

	auto *latHB = new Lattice<vSu2<Simd_f>,EvenOdd>(16, 16);
	latHB->SetRand();

	auto *latMP = new Lattice<vSu2<Simd_f>,EvenOdd>(*latHB);

	Su2Action::Action<vSu2<Simd_f>,EvenOdd> wActHB(*latHB, 2.05, 0.0);
	Su2Action::Action<vSu2<Simd_f>,EvenOdd> wActMP(*latMP, 2.05, 0.0);


	LogOut("Tuning action...\n");
	Su2Tune::Tune(wActHB);
	Su2Tune::Tune(wActMP);
	LogOut("Tuning completed\n");

	auto	plqHB = wActHB.allPts();
	auto	plqMP = wActMP.allPts();
	LogOut("\nAction value: HB: %le\tMP: %le\n\n", plqHB/(6.*latHB->Volume()), plqMP/(6.*latMP->Volume()));

	double	hbAvg = 0., hbErr = 0.;
	double	mpAvg = 0., mpErr = 0.;

	auto	sHB = Su2Action::HeatBath  <vSu2<Simd_f>,EvenOdd>(wActHB);
	auto	sOR = Su2Action::OverRelax <vSu2<Simd_f>,EvenOdd>(wActHB);
	auto	sMP = Su2Action::Metropolis<vSu2<Simd_f>,EvenOdd>(wActMP);
	auto	sPH = Su2Action::Plaquette <vSu2<Simd_f>,EvenOdd>(*latHB);
	auto	sPM = Su2Action::Plaquette <vSu2<Simd_f>,EvenOdd>(*latMP);

	LogOut("\nTuning HeatBath...\n");
	Su2Tune::Tune(sHB);
	LogOut("Tuning OverRelax...\n");
	Su2Tune::Tune(sOR);
	LogOut("Tuning Metropolis...\n");
	Su2Tune::Tune(sMP);
	LogOut("Tuning Plaquette...\n");
	Su2Tune::Tune(sPH);
	Su2Tune::Tune(sPM);
	LogOut("Tuning completed\n");

	plqHB = wActHB.allPts();
	plqMP = wActMP.allPts();
	LogOut("\nAction value: HB: %le\tMP: %le\n\n", plqHB/(6.*latHB->Volume()), plqMP/(6.*latMP->Volume()));

	LogOut("\nTermalizing...\n\n");
	start = std::chrono::high_resolution_clock::now();

	for	(int i = 0; i<nTerm; i++) {
		sHB(nOvHB);
		sOR();
		sMP(nTries);
	}

	plqHB = sPH();
	plqMP = sPM();

	LogOut("\nTermalization completed\n");
	LogOut("\nPlaquette value: HB: %le\tMP: %le\n\n", plqHB, plqMP);

	for	(int i = 0; i<nIters; i++) {
		sHB(nOvHB);
		sOR();
		sMP(nTries);

		plqHB = sPH();
		plqMP = sPM();

		hbAvg += plqHB;
		hbErr += plqHB*plqHB;

		mpAvg += plqMP;
		mpErr += plqMP*plqMP;

		LogOut("%04d\t%.7e\t%.7e\n", i, plqHB, plqMP);
	}

	stop  = std::chrono::high_resolution_clock::now();

	hbAvg /= ((double) nIters);
	hbErr /= ((double) nIters);
	mpAvg /= ((double) nIters);
	mpErr /= ((double) nIters);
	stdElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	LogOut("Elapsed time %.3lf ms\n", ((double) stdElapsed.count())/1e6);
	LogOut("Final Plaquette HB: %le +/- %le\n", hbAvg, std::sqrt((hbErr - hbAvg*hbAvg)/((double) (nIters - 1))));
	LogOut("Final Plaquette MP: %le +/- %le\n", mpAvg, std::sqrt((mpErr - mpAvg*mpAvg)/((double) (nIters - 1))));

	delete latHB;
	delete latMP;

	LogOut("\n\nTest improved action\n");

	auto *iLatHB = new Lattice<vSu2<Simd_f>,Colored>(16, 16);
	iLatHB->SetRand();

	auto *iLatMP = new Lattice<vSu2<Simd_f>,Colored>(*iLatHB);

	Su2Action::Action<vSu2<Simd_f>,Colored> iActHB(*iLatHB, 2.05, 0.0);
	Su2Action::Action<vSu2<Simd_f>,Colored> iActMP(*iLatMP, 2.05, 0.0);

	LogOut("Tuning action...\n");
	Su2Tune::Tune(iActHB);
	Su2Tune::Tune(iActMP);
	LogOut("Tuning completed\n");

	plqHB = iActHB.allPts();
	plqMP = iActMP.allPts();
	LogOut("\nAction value: HB: %le\tMP: %le\n\n", plqHB/(6.*iLatHB->Volume()), plqMP/(6.*iLatMP->Volume()));

	hbAvg = 0., hbErr = 0.;
	mpAvg = 0., mpErr = 0.;

	auto	iHB = Su2Action::HeatBath  <vSu2<Simd_f>,Colored>(iActHB);
	auto	iOR = Su2Action::OverRelax <vSu2<Simd_f>,Colored>(iActHB);
	auto	iMP = Su2Action::Metropolis<vSu2<Simd_f>,Colored>(iActMP);
	auto	iPH = Su2Action::Plaquette <vSu2<Simd_f>,Colored>(*iLatHB);
	auto	iPM = Su2Action::Plaquette <vSu2<Simd_f>,Colored>(*iLatMP);

	LogOut("\nTuning HeatBath...\n");
	Su2Tune::Tune(iHB);
	LogOut("Tuning OverRelax...\n");
	Su2Tune::Tune(iOR);
	LogOut("Tuning Metropolis...\n");
	Su2Tune::Tune(iMP);
	LogOut("Tuning Plaquette...\n");
	Su2Tune::Tune(iPH);
	Su2Tune::Tune(iPM);
	LogOut("Tuning completed\n");

	plqHB = iActHB.allPts();
	plqMP = iActMP.allPts();
	LogOut("\nAction value: HB: %le\tMP: %le\n\n", plqHB/(6.*iLatHB->Volume()), plqMP/(6.*iLatMP->Volume()));

	LogOut("\nTermalizing...\n\n");
	start = std::chrono::high_resolution_clock::now();

	for	(int i = 0; i<nTerm; i++) {
		iHB(nOvHB);
		iOR();
		iMP(nTries);
	}

	plqHB = iPH();
	plqMP = iPM();

	LogOut("\nTermalization completed\n");
	LogOut("\nPlaquette value: HB: %le\tMP: %le\n\n", plqHB, plqMP);

	for	(int i = 0; i<nIters; i++) {
		iHB(nOvHB);
		iOR();
		iMP(nTries);

		plqHB = iPH();
		plqMP = iPM();

		hbAvg += plqHB;
		hbErr += plqHB*plqHB;

		mpAvg += plqMP;
		mpErr += plqMP*plqMP;

		LogOut("%04d\t%.7e\t%.7e\n", i, plqHB, plqMP);
	}

	stop  = std::chrono::high_resolution_clock::now();

	hbAvg /= ((double) nIters);
	hbErr /= ((double) nIters);
	mpAvg /= ((double) nIters);
	mpErr /= ((double) nIters);
	stdElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	LogOut("Elapsed time %.3lf ms\n", ((double) stdElapsed.count())/1e6);
	LogOut("Final Plaquette HB: %le +/- %le\n", hbAvg, std::sqrt((hbErr - hbAvg*hbAvg)/((double) (nIters - 1))));
	LogOut("Final Plaquette MP: %le +/- %le\n", mpAvg, std::sqrt((mpErr - mpAvg*mpAvg)/((double) (nIters - 1))));

	delete	iLatHB;
	delete	iLatMP;

	endSu2  ();
}
