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

using namespace Simd;

int	main (int argc, char *argv[]) {
	initSu2 (argc, argv);
	Lattice<Su2<double>> *myLat = new Lattice<Su2<double>>(8, 8, PrecDouble);
	myLat->SetRand();

	for (int i; i<myLat->Volume(); i+=4) {
		double *ptr = myLat->RawData();
		double test = ptr[i]*ptr[i] + ptr[i+1]*ptr[i+1] + ptr[i+2]*ptr[i+2] + ptr[i+3]*ptr[i+3];
		printf("%05d\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t\t%.3lf\n", i, ptr[i], ptr[i+1], ptr[i+2], ptr[i+3], test);
	}

	Su2Action::Action<Su2<double>> wAction(*myLat, 2.0, 0.0, false);

	auto	sAct = wAction.allPts();

	printf("Plaquette value: %le\n", sAct/myLat->Volume());

	delete myLat;

	Lattice<vSu2<Simd_f>> *myVLat = new Lattice<vSu2<Simd_f>>(16, 16, PrecSingle);
	myVLat->SetRand();

	for (int i; i<myVLat->oVol(); i++) {
		for (int mu=0; mu<4; mu++) {
			auto test = myVLat->Data(i,mu);
			printf("%05d %d\t", i, mu);
			printsVar((test.Norm()).raw()," ");
		}
	}

	Su2Action::Action<vSu2<Simd_f>> wVAction(*myVLat, 2.0, 0.0, false);

	auto sVAct = wVAction.allPts();

	printf("Plaquette value: %le\n", sVAct/myVLat->Volume());

	delete myVLat;

	endSu2  ();
}
