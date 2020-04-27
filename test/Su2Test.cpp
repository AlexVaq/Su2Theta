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

#define	Ordering Colored // EvenOdd
#define Float    double  // float
#define vSmd     Simd_d  // Simd_f
using namespace Simd;

int	main (int argc, char *argv[]) {
	initSu2 (argc, argv);

	Lattice<Su2<Float>,Ordering> *myLat = new Lattice<Su2<Float>,Ordering>(16, 16);
	myLat->SetRand();

	Su2Action::Action<Su2<Float>,Ordering> wAction(*myLat, 2.0, 0.0);
	Su2Tune::Tune(wAction);
	auto	sAct = wAction.allPts();
	printf("Plaquette value: %le\n", sAct/(6.*myLat->Volume()));

	Su2Action::Plaquette<Su2<Float>,Ordering> wPlq(*myLat);
	Su2Tune::Tune(wPlq);
	printf("Plaquette value: %le\n", wPlq()/(6.*myLat->Volume()));

	printf("Scalar coordinate test: index to coordinate to index\n");
	for (size_t i=0; i<myLat->oVol(); i++) {
		Coord<Ordering> pt (i, myLat->vLength());

		//printf("Index %05zu --> %05zu Col %02d (%02d %02d %02d %02d)\n", i, pt.Index(myLat->vLength()), PointParity(pt,  *myLat), pt.x[0], pt.x[1], pt.x[2], pt.x[3]);
		if (i != pt.Index(myLat->vLength()))
			printf("\n\nVAYA CAGADA %zu %zu\n\n", i, pt.Index(myLat->vLength()));
		fflush(stdout);
	}
	printf("Done\n");

	printf("Scalar coordinate test: coordinate to index to coordinate\n");
	for (size_t it=0; it<myLat->vLength()[3]; it++) {
	  for (size_t iz=0; iz<myLat->vLength()[2]; iz++) {
	    for (size_t iy=0; iy<myLat->vLength()[1]; iy++) {
	      for (size_t ix=0; ix<myLat->vLength()[0]; ix++) {
		Coord<Ordering> pt (ix, iy, iz, it);
		auto idx = pt.Index(myLat->vLength());
		Coord<Ordering> nt (idx, myLat->vLength());
		auto clp = PointParity(pt, *myLat);
		auto cln = PointParity(nt, *myLat);

		if (pt != nt || clp != cln)
			printf("\n\nVAYA CAGADA (%d %d %d %d) vs (%d %d %d %d)\n\n", pt.x[0], pt.x[1], pt.x[2], pt.x[3], nt.x[0], nt.x[1], nt.x[2], nt.x[3]);
		fflush(stdout);
	      }
	    }
	  }
	}
	printf("Done\n");
	delete myLat;

	Lattice<vSu2<vSmd>,Ordering> *myVLat = new Lattice<vSu2<vSmd>,Ordering>(16, 16);
	myVLat->SetRand();

	Su2Action::Action<vSu2<vSmd>,Ordering> wVAction(*myVLat, 2.0, 0.0);
	Su2Tune::Tune(wVAction);

	Su2Action::Plaquette<vSu2<vSmd>,Ordering> wvPlq(*myVLat);
	Su2Tune::Tune(wvPlq);
	printf("Plaquette value: %le\n", wvPlq()/(6.*myVLat->oVol()));

	printf("Vector coordinate test: index to coordinate to index\n");
	for (size_t i=0; i<myVLat->oVol(); i++) {
		Coord<Ordering> pt (i, myVLat->vLength());

		//printf("Index %05zu --> %05zu Cold %02d (%02d %02d %02d %02d)\n", i, pt.Index(myVLat->vLength()), PointParity(pt,  *myVLat), pt.x[0], pt.x[1], pt.x[2], pt.x[3]);
		if (i != pt.Index(myVLat->vLength()))
			printf("\n\nVAYA CAGADA %zu %zu\n\n", i, pt.Index(myVLat->vLength()));
		fflush(stdout);
	}
	printf("Done\n");

	printf("Vector coordinate test: coordinate to index to coordinate\n");
	for (size_t it=0; it<myVLat->vLength()[3]; it++) {
	  for (size_t iz=0; iz<myVLat->vLength()[2]; iz++) {
	    for (size_t iy=0; iy<myVLat->vLength()[1]; iy++) {
	      for (size_t ix=0; ix<myVLat->vLength()[0]; ix++) {
		Coord<Ordering> pt (ix, iy, iz, it);
		auto idx = pt.Index(myVLat->vLength());
		Coord<Ordering> nt (idx, myVLat->vLength());
		auto clp = PointParity(pt, *myVLat);
		auto cln = PointParity(nt, *myVLat);

		if (pt != nt || clp != cln)
			printf("\n\nVAYA CAGADA (%d %d %d %d) vs (%d %d %d %d)\n\n", pt.x[0], pt.x[1], pt.x[2], pt.x[3], nt.x[0], nt.x[1], nt.x[2], nt.x[3]);
		//printf("\n\n(%d %d %d %d | %d) vs (%d %d %d %d | %d)\n\n", pt.x[0], pt.x[1], pt.x[2], pt.x[3], clp, nt.x[0], nt.x[1], nt.x[2], nt.x[3], cln);
		fflush(stdout);
	      }
	    }
	  }
	}
	printf("Done\n");
	auto sVAct = wVAction.allPts();

	printf("Vector Plaquette value: %le\n", sVAct/(6.*myVLat->Volume()));

	Unfolder<vSu2<vSmd>,Ordering>	munger(*myVLat);

	Lattice<Su2<Float>,Ordering> myULat = munger();

	Su2Action::Action<Su2<Float>,Ordering> wUAction(myULat, 2.0, 0.0);
	auto sUAct = wUAction.allPts();

	printf("Unfold Plaquette value: %le\n", sUAct/(6.*myULat.Volume()));
/*
	double	acc = 0.;
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
	Lattice<vSu2<vSmd>> *myDLat = new Lattice<vSu2<vSmd>>(16, 16);
	myDLat->SetRand();

	Su2Action::Action<vSu2<vSmd>> wDAction(*myDLat, 2.0, 0.0, false);
	sVAct = wDAction.allPts();
	printf("Plaquette value: %le\n", sVAct/(6.*myDLat->Volume()));

	Unfolder	dMunger(*myDLat);

	Lattice<Su2<Float>> mySLat = dMunger();

	Su2Action::Action<Su2<Float>> wSAction(mySLat, 2.0, 0.0, false);
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
