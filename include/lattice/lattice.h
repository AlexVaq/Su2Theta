#include <enumFields.h>
#include "utils/memAlloc.h"
#include "simd/simd.h"

using namespace Su2Enum;

template<class T>
class	Lattice	{

	private:

	int		Ls, Lt;
	int		sVol, tVol;
	Device		gDev;

	double		beta;

	void		*hData;
	void		*hMoms;

	size_t		dSize;

	public:

		Lattice<T> (int Ls, int Lt, Precision prec) : Ls(Ls), Lt(Lt), sVol(Ls*Ls*Ls), tVol(Ls*Ls*Ls*Lt), dSize(sizeof(T)) {

		size_t  dataSize = tVol*dSize;

		#ifndef USE_GPU
		alignAlloc((void **) &hData, Simd::sAlign, dataSize);
		#else
		if (gDev == DeviceGpuManaged) {
			cudaMallocManaged (&gData, dataSize);
			hData = gData;
		} else {
			alignAlloc(&hData, Simd::sAlign, dataSize);

			if (gDev == DeviceGpu)
				cudaMalloc (&gData, dataSize);
			CudaCheckErrors();
		}
		#endif

//		hMoms = ALLOCATE MEMORY;
	}

		~Lattice<T> () {
		#ifndef USE_GPU
		if (hData != nullptr)
			trackFree(hData);
		#else
		if (gDev == DeviceGpuManaged) {
			if (gData != nullptr)
				cudaFree(gData);
		} else {
			if (hData != nullptr)
				trackFree(hData);

			if (gDev == DeviceGpu && gData != nullptr)
				cudaFree(gData);
		}

		CudaCheckErrors();
		#endif
//		if (hMoms != nullptr)
//			trackFree(hMoms);
	}

	void	SetRand	();

	int	SLength	() { return Ls; }
	int	TLength	() { return Lt; }

	int	Volume	() { return tVol; }
	int	SVol	() { return sVol; }
};
