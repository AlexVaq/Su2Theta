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

	Su2Action::Action<Su2<double>> wAction(*myLat, 2.0, 0.0, false);

	auto	sAct = wAction.allPts();

	printf("Plaquette value: %le\n", sAct/myLat->Volume());

	delete myLat;

	Lattice<vSu2<Simd_f>> *myVLat = new Lattice<vSu2<Simd_f>>(16, 16, PrecSingle);
	myVLat->SetRand();

	Su2Action::Action<vSu2<Simd_f>> wVAction(*myVLat, 2.0, 0.0, false);

	auto sVAct = wVAction.allPts();

	printf("Plaquette value: %le\n", sVAct/myVLat->Volume());

	delete myVLat;

	endSu2  ();
}
