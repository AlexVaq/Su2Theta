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
#include "utils/unfolder.h"

using namespace Simd;

int	main (int argc, char *argv[]) {
	initSu2 (argc, argv);

	Lattice<Su2<float>,EvenOdd> *myLat = new Lattice<Su2<float>,EvenOdd>(16, 16);
	myLat->SetRand();

	Su2Action::Action<Su2<float>,EvenOdd> wAction(*myLat, 2.0, 0.0);
	Su2Tune::Tune(wAction);
	auto	sAct = wAction.allPts();
	printf("Plaquette value: %le\n", sAct/(6.*myLat->Volume()));

	Su2Action::Plaquette<Su2<float>,EvenOdd> wPlq(*myLat);
	Su2Tune::Tune(wPlq);
	printf("Plaquette value: %le\n", wPlq()/(6.*myLat->Volume()));

	printf("Coordinate test\n");
	for (size_t i=0; i<myLat->oVol(); i++) {
		Coord<EvenOdd> pt (i, myLat->vLength());

		//printf("Index %05zu --> %05zu Col %02d (%02d %02d %02d %02d)\n", i, pt.Index(myLat->vLength()), PointParity(pt,  *myLat), pt.x[0], pt.x[1], pt.x[2], pt.x[3]);
		if (i != pt.Index(myLat->vLength()))
			printf("\n\nVAYA CAGADA %zu %zu\n\n", i, pt.Index(myLat->vLength()));
		//fflush(stdout);
	}
	printf("Done\n");

	delete myLat;

	Lattice<vSu2<Simd_f>,EvenOdd> *myVLat = new Lattice<vSu2<Simd_f>,EvenOdd>(16, 16);
	myVLat->SetRand();

	Su2Action::Action<vSu2<Simd_f>,EvenOdd> wVAction(*myVLat, 2.0, 0.0);
	Su2Tune::Tune(wVAction);

	Su2Action::Plaquette<vSu2<Simd_f>,EvenOdd> wvPlq(*myVLat);
	Su2Tune::Tune(wvPlq);
	printf("Plaquette value: %le\n", wvPlq()/(6.*myVLat->oVol()));

	printf("Coordinate test\n");
	for (size_t i=0; i<myVLat->oVol(); i++) {
		Coord<EvenOdd> pt (i, myVLat->vLength());

		//printf("Index %zu --> %zu (%02d %02d %02d %02d)\n", i, pt.Index(myVLat->vLength()), pt.x[0], pt.x[1], pt.x[2], pt.x[3]);
		if (i != pt.Index(myVLat->vLength()))
			printf("\n\nVAYA CAGADA %zu %zu\n\n", i, pt.Index(myVLat->vLength()));
		//fflush(stdout);
	}
	printf("Done\n");

	auto sVAct = wVAction.allPts();

	printf("Vector Plaquette value: %le\n", sVAct/(6.*myVLat->Volume()));

	Unfolder<vSu2<Simd_f>,EvenOdd>	munger(*myVLat);

	Lattice<Su2<float>,EvenOdd> myULat = munger();

	Su2Action::Action<Su2<float>,EvenOdd> wUAction(myULat, 2.0, 0.0);
	auto sUAct = wUAction.allPts();

	printf("Unfold Plaquette value: %le\n", sUAct/(6.*myULat.Volume()));
/*
	float	acc = 0.;
	for (size_t i=0; i<myVLat->oVol(); i++) {
		for (size_t j=0; j<4; j++) {
			auto stp = wVAction(i,j);
			auto lnk = myVLat->Data(i,j);
			acc += (lnk*stp).SuperTrace();
		}
	}

	printf("Acc %le\n", acc/(72.*myVLat->Volume()));
*/
	delete myVLat;
/*
	Lattice<vSu2<Simd_f>> *myDLat = new Lattice<vSu2<Simd_f>>(16, 16);
	myDLat->SetRand();

	Su2Action::Action<vSu2<Simd_f>> wDAction(*myDLat, 2.0, 0.0, false);
	sVAct = wDAction.allPts();
	printf("Plaquette value: %le\n", sVAct/(6.*myDLat->Volume()));

	Unfolder	dMunger(*myDLat);

	Lattice<Su2<float>> mySLat = dMunger();

	Su2Action::Action<Su2<float>> wSAction(mySLat, 2.0, 0.0, false);
	sUAct = wSAction.allPts();
	printf("Plaquette value: %le\n", sUAct/(6.*mySLat.Volume()));

	Su2Action::OverRelax(wDAction);
	Su2Action::OverRelax(wSAction);

	sVAct = wDAction.allPts();
	printf("Plaquette value: %le\n", sVAct/(6.*myDLat->Volume()));
	sUAct = wSAction.allPts();
	printf("Plaquette value: %le\n", sUAct/(6.*mySLat.Volume()));

	delete	myDLat;
*/
	endSu2  ();
}
