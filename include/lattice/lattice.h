#ifndef	__LATTICECLASS
	#define	__LATTICECLASS

	#include <enumFields.h>
	#include "utils/memAlloc.h"
	#include "simd/simd.h"

	using namespace Su2Enum;

	struct	Point {
		size_t		idx;
		Permutation	perm;

		Point (size_t idx, Permutation perm) : idx(idx), perm(perm) {}
		Point (size_t idx) : idx(idx), perm(NoPermutation) {}
		Point () : idx(0), perm(NoPermutation) {}
	};

	template<class T>
	class	Lattice	{

		private:

		size_t		Ls,  Lt;
		size_t		vLx, vLy, vLz, vLt;
		size_t		sVol,  tVol;
		size_t		vsVol, vtVol;
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

			vLx   = Ls/T::xWide;
			vLy   = Ls/T::yWide;
			vLz   = Ls/T::zWide;
			vLt   = Ls/T::tWide;
			vsVol = vLx*vLy*vLz;
			vtVol = vsVol*vLt;
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

		size_t		SLength	() { return Ls;  }
		size_t		vxLength() { return vLx; }
		size_t		vyLength() { return vLy; }
		size_t		vzLength() { return vLz; }
		size_t		vtLength() { return vLt; }
		size_t		TLength	() { return Lt;  }

		size_t		Volume	() { return tVol; }
		size_t		SVol	() { return sVol; }

		size_t		oVol	() { return vtVol; }
		size_t		oSVol	() { return vsVol; }

		Point		nextPoint (Point pt, int mu) {
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
				case	0:	// X+ direction
				{
					auto x0 = pt.idx%vLx;

					if (x0 == vLx - 1) {
						pt.idx -= vLx - 1;
						pt.perm ^= PermutationX;
					} else
						pt.idx++;
				}
				break;

				case	1:	// Y+ direction
				{
					auto S  = vLx*vLy;
					auto x1 = (pt.idx/vLx)%vLy;

					if (x1 == vLy - 1) {
						pt.idx += vLx - S;
						pt.perm ^= PermutationY;
					} else
						pt.idx += vLx;
				}
				break;

				case	2:	// Z+ direction
				{
					auto S  = vLx*vLy;
					auto x2 = (pt.idx/S)%vLz;

					if (x2 == vLz - 1) {
						pt.idx += S - vsVol;
						pt.perm ^= PermutationZ;
					} else
						pt.idx += S;
				}
				break;

				case	3:	// T+ direction
				{
					auto x3 = pt.idx/vsVol;

					if (x3 == vLt - 1) {
						pt.idx += vsVol - vtVol;
						pt.perm ^= PermutationT;
					} else
						pt.idx += vsVol;
				}
				break;

				case	4:	// X- direction
				{
					auto x0 = pt.idx%vLx;

					if (x0 == 0) {
						pt.idx += (vLx - 1);
						pt.perm ^= PermutationX;
					} else
						pt.idx--;
				}
				break;

				case	5:	// Y- direction
				{
					auto S  = vLx*vLy;
					auto x1 = (pt.idx/vLx)%vLy;

					if (x1 == 0) {
						pt.idx += (S - vLx);
						pt.perm ^= PermutationY;
					} else
						pt.idx -= vLx;
				}
				break;

				case	6:	// Z- direction
				{
					auto S  = vLx*vLy;
					auto x2 = (pt.idx/S)%vLz;

					if (x2 == 0) {
						pt.idx += (vsVol - S);
						pt.perm ^= PermutationZ;
					} else
						pt.idx -= S;
				}
				break;

				case	7:	// T- direction
				{
					auto x3 = pt.idx/vsVol;

					if (x3 == 0) {
						pt.idx += (vtVol - vsVol);
						pt.perm ^= PermutationT;
					} else
						pt.idx -= vsVol;
				}
				break;
			}

			return	pt;
		}

		typename T::data*	RawData	() 	 { return static_cast<typename T::data*>      (hData); }
		const typename T::data*	RawData	() const { return static_cast<const typename T::data*>(hData); }

		T		Data	(Point pt, size_t mu)	      { T tmp =	static_cast<T*>(hData)[mu*tVol+pt.idx];
									if (pt.perm & PermutationX) tmp = tmp.xPermute();
									if (pt.perm & PermutationY) tmp = tmp.yPermute();
									if (pt.perm & PermutationZ) tmp = tmp.zPermute();
									if (pt.perm & PermutationT) tmp = tmp.tPermute();
									return	tmp; }
		void		Insert	(T tmp, size_t idx, size_t mu)  { static_cast<T*>(hData)[mu*tVol+idx] = tmp; }
	};

	template<class T>
	ParityType	PointParity (size_t idx, Lattice<T> &myLat) {
		auto vLx   = myLat.vxLength();
		auto vLy   = myLat.vyLength();
		auto vLz   = myLat.vzLength();
		auto S     = vLx*vLy;
		auto vsVol = myLat.oSVol();

		auto x0 = idx%vLx;
		auto x1 = (idx/vLx)%vLy;
		auto x2 = (idx/S)%vLz;
		auto x3 = idx/vsVol;

		return	((ParityType)((x0 + x1 + x2 + x3) & 1));
	}

	template<class T>
	ParityType	PointParity (Point pt, Lattice<T> &myLat) {
		return	PointParity(pt.idx, myLat);
	}
#endif
