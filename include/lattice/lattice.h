#ifndef	__LATTICECLASS
	#define	__LATTICECLASS

	#include <enumFields.h>
	#include "utils/memAlloc.h"
	#include "simd/simd.h"

	using namespace Su2Enum;

	template<class T>
	class	Lattice	{

		private:

		size_t		Ls, Lt;
		size_t		sVol, tVol;
		Device		gDev;

		double		beta;

		void		*hData;
		void		*hMoms;

		size_t		dSize;

		public:

			Lattice<T> (int Ls, int Lt, Precision prec) : Ls(Ls), Lt(Lt), sVol(Ls*Ls*Ls), tVol(Ls*Ls*Ls*Lt), dSize(sizeof(T)) {

			size_t  dataSize = 4*tVol*dSize;

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

//			hMoms = ALLOCATE MEMORY;
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
//			if (hMoms != nullptr)
//				trackFree(hMoms);
		}

		void	SetRand	() {
			// GET THE PROFILER
			#pragma omp parallel for schedule(static)
			for (int i=0; i<tVol*4; i++)
				static_cast<T*>(hData)[i].SetRandom();
		}

		size_t		SLength	() { return Ls; }
		size_t		TLength	() { return Lt; }

		size_t		Volume	() { return tVol; }
		size_t		SVol	() { return sVol; }

		size_t		oVol	() { return tVol/T::sWide; }

		size_t		nextPoint (size_t idx, int mu) {
			/*	Even-Odd code	*/
#if 0
			size_t	x[4];
			size_t	hVol   = tVol >> 1;
			size_t	parity = idx > hVol ? 1 : 0;

			size_t	tmp1   = idx  / (Ls >> 1);
			size_t	tmp2   = tmp1 / Ls;
			size_t  x[1]   = tmp1 - tmp2 * Ls;
			size_t  x[3]   = tmp2 / Ls;
			size_t  x[2]   = tmp2 - x[3] * Ls;
			size_t  x[0]   = 2*idx + ((x[1] + x[2] + x[3] + parity) & 1) - tmp1 * Ls;

			if (mu != 3) {
				if (x[mu] == Ls - 1)
					x[mu] -= Ls - 1;
				else
					x[mu]++;
			} else {
				if (x[mu] == Lt - 1)
					x[mu] -= Lt - 1;
				else
					x[mu]++;
			}

			return (x[0] >> 1) + (Ls >> 1) * x[1] + (Ls >> 1)*Ls * x[2] + (sVol >> 1) * x[3] + ((x[0] + x[1] + x[2] + x[3]) & 1) * hVol;
#endif
			/*	Lexicographic code	*/
			switch (mu) {
				case	0:	// X direction
				{
					auto x0 = idx%Ls;

					if (x0 == Ls - 1)
						return	idx - Ls + 1;
					else
						return	idx + 1;
				}
				break;

				case	1:	// Y direction
				{
					auto S  = Ls*Ls;
					auto x1 = (idx/Ls)%Ls;

					if (x1 == Ls - 1)
						return	idx + Ls - S;
					else
						return	idx + Ls;
				}
				break;

				case	2:	// Z direction
				{
					auto x2 = (idx/(Ls*Ls))%Ls;

					if (x2 == Ls - 1)
						return	idx + Ls*Ls - sVol;
					else
						return	idx + Ls*Ls;
				}
				break;

				case	3:	// T direction
				{
					auto x3 = idx/sVol;

					if (x3 == Lt - 1)
						return	idx + sVol - tVol;
					else
						return	idx + sVol;
				}
				break;
			}
		}

		typename T::data*	RawData	() 	 { return static_cast<typename T::data*>      (hData); }
		const typename T::data*	RawData	() const { return static_cast<const typename T::data*>(hData); }

		T		Data	(size_t idx, size_t mu)	      { return	static_cast<T*>      (hData)[mu*tVol+idx]; }
//		template<const bool dagger=false>
//		const T&	Data	(size_t idx, size_t mu) const { return dagger ? !(static_cast<const T*>(hData)[mu*tVol+idx]) : static_cast<const T*>(hData)[mu*tVol+idx]; }

		T		nextData(size_t idx, size_t mu, int dir) {
			size_t	rIdx;
			/*	Lexicographic code	*/
			switch (dir) {
				case	0:	// X direction
				{
					auto x0 = idx%Ls;

					if (x0 == Ls - 1)
						return	Data(idx + 1 - Ls, mu).xPermute();
					else
						return	Data(idx + 1, mu);
				}
				break;

				case	1:	// Y direction
				{
					auto S  = Ls*Ls;
					auto x1 = (idx/Ls)%Ls;

					if (x1 == Ls - 1)
						return	Data(idx + Ls - S, mu).yPermute();
					else
						return	Data(idx + Ls, mu);
				}
				break;

				case	2:	// Z direction
				{
					auto x2 = (idx/(Ls*Ls))%Ls;

					if (x2 == Ls - 1)
						return	Data(idx + Ls*Ls - sVol, mu).zPermute();
					else
						return	Data(idx + Ls*Ls, mu);
				}
				break;

				case	3:	// T direction
				{
					auto x3 = idx/sVol;

					if (x3 == Lt - 1)
						return	Data(idx + sVol - tVol, mu).tPermute();
					else
						return	Data(idx + sVol, mu);
				}
				break;
			}
		}
	};
#endif
