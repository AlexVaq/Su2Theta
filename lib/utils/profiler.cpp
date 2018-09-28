#include <string>
#include <vector>
#include <map>
#include "utils/utils.h"
#include "enumFields.h"
#include "comms/comms.h"

using namespace	Su2Enum;

namespace Su2Prof {

	static	std::map<ProfType,Profiler>			profs;
	static	std::chrono::high_resolution_clock::time_point	stPt;

	static	double tTime = 0.;

	double	Profiler::printStats	() {
		double	aTime = 0.;

		for (auto data = prof.cbegin(); data != prof.cend(); data++)
	        {
			std::string	name   = data->first;
		        FlopCounter	fCount = data->second;

			aTime += fCount.DTime();

			LogMsg (VerbSilent, "\tFunction %-30s\tGFlops %.3lf\tGBytes %.3lf\tTotal time %.2lfs (%.2lf\%)", name.c_str(), fCount.GFlops(), fCount.GBytes(), fCount.DTime(), 100.*fCount.DTime()/tTime);
        	}

		return	aTime;
	}

	void	initProfilers() {

		stPt = std::chrono::high_resolution_clock::now();

		Profiler	genProfiler("Genconf");
		profs.insert(std::make_pair(ProfGen, genProfiler));

		Profiler	actProfiler("Action");
		profs.insert(std::make_pair(ProfAction, actProfiler));

		Profiler	plaqProfiler("Plaquette");
		profs.insert(std::make_pair(ProfPlaq, plaqProfiler));

		Profiler	hbProfiler("HeatBath");
		profs.insert(std::make_pair(ProfHB, hbProfiler));

		Profiler	ovrProfiler("OverRelax");
		profs.insert(std::make_pair(ProfOvR, ovrProfiler));

		Profiler	mpProfiler("Metro");
		profs.insert(std::make_pair(ProfMetro, mpProfiler));

		Profiler	qProfiler("QCharge");	// Falta escribir
		profs.insert(std::make_pair(ProfQCharge, qProfiler));

		Profiler	tunerProfiler("Tuner");
		profs.insert(std::make_pair(ProfTuner, tunerProfiler));

		Profiler	folderProfiler("Folder");
		profs.insert(std::make_pair(ProfFold, folderProfiler));

		Profiler	hdf5Profiler("Hdf5 I/O");	// Falta escribir
		profs.insert(std::make_pair(ProfHdf5, hdf5Profiler));
	}

	void	printProfStats() {
		if (Su2Comms::rank != 0)
			return;

		double	aTime = 0.;

		tTime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - stPt).count()*1e-6;

		for (auto &data : profs) {
			auto &cProf = data.second;
			LogMsg(VerbSilent, "Profiler %s:", cProf.name().c_str());
			auto cTime = cProf.printStats();
			aTime += cTime;
			LogMsg(VerbSilent, "Total %s: %.2lf", cProf.name().c_str(), cTime);
			LogMsg(VerbSilent, "");
		}
		LogMsg (VerbSilent, "Unaccounted time %.2lfs of %.2lfs (%.2lf\%)", tTime - aTime, tTime, 100.*(1. - aTime/tTime));
	}

	Profiler&	getProfiler(ProfType pType) {
		return	profs[pType];
	}
};
