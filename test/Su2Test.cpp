#include <memory>
#include "enumFields.h"
#include "lattice/lattice.h"
#include "su2/su2.h"
#include "utils/memAlloc.h"
#include "utils/logger.h"
#include "utils/profiler.h"
#include "utils/misc.h"

int	main (int argc, char *argv[]) {
	initSu2 (argc, argv);
	std::unique_ptr<Lattice<Su2<double>>> myLat = std::make_unique<Lattice<Su2<double>>>(8, 8, PrecDouble);
	myLat.reset();	// For the logger
	endSu2  ();
}
