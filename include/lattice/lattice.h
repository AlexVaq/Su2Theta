#ifndef	__LATTICECLASS
	#define	__LATTICECLASS

	#include <enumFields.h>
	#include "utils/utils.h"
	#include "simd/simd.h"
//#define	Lexi
	using namespace Su2Enum;

	struct	Point {
		size_t		idx;
		Permutation	perm;

		Point (size_t idx, Permutation perm) : idx(idx), perm(perm) {}
		Point (size_t idx) : idx(idx), perm(NoPermutation) {}
		Point () : idx(0), perm(NoPermutation) {}
	};

	template<const ColorKind cOrd>
	struct	Coord {
		int		x[4];
		Permutation	perm;

		Coord (int x0, int x1, int x2, int x3, Permutation perm) : x({ x0, x1, x2, x3 }), perm(perm) {}
		Coord (int x0, int x1, int x2, int x3)			 : x({ x0, x1, x2, x3 }), perm(NoPermutation) {}
		Coord (size_t idx, size_t vL[4])	 		 : perm(NoPermutation) {

			switch (cOrd) {
				case EvenOdd: {
#ifdef	Lexi
					auto tmp = idx/vL[0];
	     				x[0] = idx%vL[0];
					x[1] = tmp%vL[1];
					tmp /= vL[1];
					x[2] = tmp%vL[2];
					x[3] = tmp/vL[2];
#else
					auto xR = (vL[0]>>1);
					auto hW = xR*vL[1]*vL[2]*vL[3];

					ParityType pBase = (idx >= hW) ? ParityOdd : ParityEven;

					auto ix = idx - hW*pBase;
					auto za = ix/xR;
					auto zb = za/vL[1];

					x[1] = (za - zb*vL[1]);
					x[3] = (zb/vL[2]);
					x[2] = (zb - x[3]*vL[2]);
					auto xO = (x[1] + x[2] + x[3] + pBase)&1;
					x[0] = 2*ix + xO - za*vL[0];
#endif
				}
				break;

				case Colored: {
#ifdef	Lexi
					auto tmp = idx/vL[0];
	     				x[0] = idx%vL[0];
					x[1] = tmp%vL[1];
					tmp /= vL[1];
					x[2] = tmp%vL[2];
					x[3] = tmp/vL[2];
#else
					auto xR = (vL[0]>>2);
					auto yR = (vL[1]>>1);
					auto zR = (vL[2]>>1);
					auto hW = (vL[3]>>1)*zR*yR*xR;

					ParityColor pCol  = (ParityColor) (idx/hW);
					Parity2Type p2    = (Parity2Type) (pCol&1);
					ParityType  pBase = (pCol >= 16) ? ParityOdd : ParityEven;

					auto ix = (idx - hW*pCol);

					auto za = ix/xR;
					auto zb = za/yR;

					x[3] = ((zb/zR)<<1) + ((pCol&8)>>3);
					x[2] = (zb<<1) - (x[3]>>1)*vL[2] + ((pCol&4)>>2);
					x[1] = (za<<1) - zb*vL[2] + ((pCol&2)>>1);
					x[0] = (ix<<2) - za*vL[0] + (((p2<<1) + (x[1]&2) + (x[2]&2) + (x[3]&2))&2) + ((pBase + (x[1]&1) + (x[2]&1) + (x[3]&1))&1);
#endif
				}
				break;
			}
		}
		Coord () : x({ 0, 0, 0, 0 }), perm(NoPermutation) {}

		inline	size_t	Index(size_t vL[4]) {

			switch (cOrd) {
				case EvenOdd: {
#ifdef	Lexi
					return	x[0] + vL[0]*(x[1] + vL[1]*(x[2] + vL[2]*x[3]));
#else
					ParityType  pBase;
					pBase = (ParityType)  ((  x[0]      +  x[1]      +  x[2]      +  x[3]     ) & 1);
					return	(x[0]>>1) + (vL[0]>>1)*(x[1] + vL[1]*(x[2] + vL[2]*(x[3] + pBase*vL[3])));
#endif
				}
				break;

				case Colored: {
#ifdef	Lexi
					return	x[0] + vL[0]*(x[1] + vL[1]*(x[2] + vL[2]*x[3]));
#else
					auto hL =    (vL[0]>>2);
					auto hS = hL*(vL[1]>>1);
					auto hV = hS*(vL[2]>>1);
					auto hW = hV*(vL[3]>>1);

					ParityColor pCol;
					ParityType  pBase;
					Parity2Type p2;

					pBase = (ParityType)  ((  x[0]      +  x[1]      +  x[2]      +  x[3]     ) & 1);
					p2    = (Parity2Type) ((((x[0] & 2) + (x[1] & 2) + (x[2] & 2) + (x[3] & 2)) & 2) >> 1);

					pCol  = (ParityColor) (p2 + ((x[1]&1)<<1) + ((x[2]&1)<<2) + ((x[3]&1)<<3) + (pBase<<4));
					// Falta xodd y x2odd	
					return	(x[0]>>2) + (x[1]>>1)*hL + (x[2]>>1)*hS + (x[3]>>1)*hV + pCol*hW;
#endif
				}
				break;
			}
		}

		inline	Point	IdxPt(size_t vL[4]) {
			return	Point(Index(vL), perm);
		}

		bool	operator==(const Coord &oPt) {
			if (oPt.x[0] != x[0])
				return false;

			if (oPt.x[1] != x[1])
				return false;

			if (oPt.x[2] != x[2])
				return false;

			if (oPt.x[3] != x[3])
				return false;

			return true;
		}

		bool	operator!=(const Coord &oPt) {
			return !((*this)==oPt);
		}
	};

	template<class T, ColorKind cOrd = EvenOdd>
	class	Lattice	: public Tunable {

		private:

		//size_t		Ls,  Lt; Inherited from Tunable
		size_t		vLx, vLy, vLz, vLt;
		size_t		vL[4];
		size_t		sVol,  tVol;
		size_t		vsVol, vtVol;
		Device		gDev;

		double		beta;

		void		*hData;
		void		*hMoms;

		size_t		dSize;
		static constexpr Permutation Perms[4] = { PermutationX, PermutationY, PermutationZ, PermutationT };

		public:

			Lattice<T,cOrd> (int Ls, int Lt) : sVol(Ls*Ls*Ls), tVol(Ls*Ls*Ls*Lt), dSize(sizeof(T)) {

			SetVolume(Ls, Lt);
			size_t  dataSize = 4*tVol*dSize/T::sWide;

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
			vLt   = Lt/T::tWide;
			vL[0] = vLx; vL[1] = vLy; vL[2] = vLz; vL[3] = vLt;
			vsVol = vLx*vLy*vLz;
			vtVol = vsVol*vLt;

			SetName  (std::string("LatticeGen"));
			SetColor (std::string(cOrd == Su2Enum::EvenOdd ? "EO" : "P32"));
			SetPrec  (sizeof(typename T::sData) == 4 ? "float" : "double");
			SetVolume(Ls, Lt);
			/*
			std::string sName = std::to_string(T::sWide) + std::string(" ") + std::to_string(Ls) + std::string("x") + std::to_string(Lt);
			if (cOrd == Su2Enum::EvenOdd)
				SetName(std::string("LatticeGen EO\t")  + sName);
			else
				SetName(std::string("LatticeGen P32\t") + sName);
			*/
		}

			Lattice<T,cOrd>	(const Lattice<T,cOrd> &myLat) {

			dSize = sizeof(T);

			size_t  dataSize = 4*myLat.tVol*dSize/T::sWide;

			#ifndef USE_GPU
			alignAlloc((void **) &hData, Simd::sAlign, dataSize);
			memcpy(hData, myLat.hData, dataSize);
			#else
			if (gDev == DeviceGpuManaged) {
				cudaMallocManaged (&gData, dataSize);
				hData = gData;
				memcpy(hData, myLat.hData, dataSize);
			} else {
				alignAlloc(&hData, Simd::sAlign, dataSize);

				if (gDev == DeviceGpu)
					cudaMalloc (&gData, dataSize);
				cudaMemcpy(gData, myLat.gData, dataSize, cudaMemcpyDeviceToDevice);
				CudaCheckErrors();
			}
			#endif

//			hMoms = ALLOCATE MEMORY;

			Ls = myLat.Ls; Lt = myLat.Lt; tVol = myLat.tVol; sVol = myLat.sVol;

			vLx   = Ls/T::xWide;
			vLy   = Ls/T::yWide;
			vLz   = Ls/T::zWide;
			vLt   = Lt/T::tWide;
			vL[0] = vLx; vL[1] = vLy; vL[2] = vLz; vL[3] = vLt;
			vsVol = vLx*vLy*vLz;
			vtVol = vsVol*vLt;

			SetName  (myLat.Name ());
			SetColor (myLat.Color());
			SetPrec  (myLat.Prec ());
			//SetVolume(myLat.SLength(), myLat.TLength());
		}

			~Lattice<T,cOrd> () {
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
			Su2Prof::Profiler &prof = Su2Prof::getProfiler(ProfGen);
			prof.start();

			#pragma omp parallel for schedule(static)
			for (int i=0; i<vtVol*4; i++)
				static_cast<T*>(hData)[i].SetRandom();

			prof.stop();
			double myGFlops = 64.0 * tVol * 1e-9;
			double myGBytes =  4.0 * vtVol * sizeof(T) / 1073741824.0;

			add(myGFlops, myGBytes); 
			prof.add(FullName(), myGFlops, myGBytes);

			LogMsg  (VerbHigh, "%s reporting %lf GFlops %lf GBytes", FullName().c_str(), prof.Prof()[FullName()].GFlops(), prof.Prof()[FullName()].GBytes());
		}

		void	SetConst	() {
			Su2Prof::Profiler &prof = Su2Prof::getProfiler(ProfGen);
			prof.start();

			//auto baseMatrix = static_cast<T*>(hData)[0];
			//baseMatrix.SetRandom();
			T baseMatrix(0.5, 0.5, 0.5, 0.5);

			#pragma omp parallel for schedule(static)
			for (int i=0; i<vtVol*4; i++)
				static_cast<T*>(hData)[i] = baseMatrix;

			prof.stop();
			double myGFlops = 0.0;
			double myGBytes = 4.0 * vtVol * sizeof(T) / 1073741824.0;

			add(myGFlops, myGBytes); 
			prof.add(FullName(), myGFlops, myGBytes);

			LogMsg  (VerbHigh, "%s reporting %lf GFlops %lf GBytes", FullName().c_str(), prof.Prof()[FullName()].GFlops(), prof.Prof()[FullName()].GBytes());
		}

		inline	size_t	SLength	() { return Ls;  }
		inline	size_t	TLength	() { return Lt;  }
		inline	size_t	vxLength() { return vLx; }
		inline	size_t	vyLength() { return vLy; }
		inline	size_t	vzLength() { return vLz; }
		inline	size_t	vtLength() { return vLt; }
		inline	size_t*	vLength()  { return vL; }

		inline	size_t	Volume	() { return tVol; }
		inline	size_t	SVol	() { return sVol; }

		inline	size_t	oVol	() { return vtVol; }
		inline	size_t	oSVol	() { return vsVol; }

		inline	Coord<cOrd>	nextPoint (Coord<cOrd> pt, int mu) {

			if (mu < 4) {
				if (pt.x[mu] == vL[mu] - 1) {
					pt.x[mu] = 0;
					pt.perm ^= Perms[mu];
				} else
					pt.x[mu]++;
			} else {
				auto dr = mu%4;
				//printf("Mu %d Dr %d x %d ", mu, dr, pt.x[dr]);
				if (pt.x[dr] == 0) {
					pt.x[dr] = vL[dr] - 1;
					pt.perm ^= Perms[dr];
				} else
					pt.x[dr]--;
				//printf("y %d\n", pt.x[dr]);
			}

			return	pt;
		}
#if 0
		inline	Point	nextPoint (Point pt, int mu) {
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
#endif
		inline	typename T::data*	RawData	() 	 { return static_cast<      typename T::data*>(hData); }
		inline	const typename T::data*	RawData	() const { return static_cast<const typename T::data*>(hData); }

		inline	T	Data	(Coord<cOrd> pt, size_t mu) { size_t idx = pt.Index(vL);
								T tmp(&(static_cast<typename T::sData*>(hData)[(mu*vtVol+idx)*sizeof(T)/sizeof(typename T::sData)]));
								if (pt.perm & PermutationX) tmp = tmp.xPermute(); 
								if (pt.perm & PermutationY) tmp = tmp.yPermute(); 
								if (pt.perm & PermutationZ) tmp = tmp.zPermute(); 
								if (pt.perm & PermutationT) tmp = tmp.tPermute(); 
								return tmp; }
		inline	T	Data	(Point pt, size_t mu) { T tmp(&(static_cast<typename T::sData*>(hData)[(mu*vtVol+pt.idx)*sizeof(T)/sizeof(typename T::sData)]));
								if (pt.perm & PermutationX) tmp = tmp.xPermute(); 
								if (pt.perm & PermutationY) tmp = tmp.yPermute(); 
								if (pt.perm & PermutationZ) tmp = tmp.zPermute(); 
								if (pt.perm & PermutationT) tmp = tmp.tPermute(); 
								return tmp; }
		inline	void	Insert 	(T &&tmp, size_t idx, size_t mu) { tmp.Stream(&(static_cast<typename T::sData*>(hData)[(mu*vtVol+idx)*sizeof(T)/sizeof(typename T::sData)])); }
		inline	void	InsertMask(typename T::Mask msk, T &&tmp, size_t idx, size_t mu) { tmp.SaveMask(msk, &(static_cast<typename T::sData*>(hData)[(mu*vtVol+idx)*sizeof(T)/sizeof(typename T::sData)])); }


		void		SaveState()	override	{
			size_t  dataSize = 4*tVol*dSize/T::sWide;

			trackAlloc((void **) &backup, dataSize);

			#ifndef USE_GPU
			memcpy(backup, hData, dataSize);
			#else
			if (gDev == DeviceGpuManaged)
				memcpy(backup, gData, dataSize);
			else
				cudaMemcpy(backup, gData, dataSize, cudaMemcpyDeviceToHost);
			#endif
		}

		void		RestoreState()	override	{
			size_t  dataSize = 4*tVol*dSize/T::sWide;

			#ifndef USE_GPU
			memcpy(hData, backup, dataSize);
			#else
			if (gDev == DeviceGpuManaged)
				memcpy(gData, backup, dataSize);
			else
				cudaMemcpy(gData, backup, dataSize, cudaMemcpyHostToDevice);
			#endif

			trackFree(backup);
		}

		inline void	Run()	override	{ SetRand(); }
		inline void	Reset() override	{ Su2Prof::getProfiler(ProfGen).reset(FullName()); }
	};

	template<class T>
	inline	auto	PointParity(Coord<EvenOdd> pt, Lattice<T,EvenOdd> &myLat) {
		return	(ParityType)  ((pt.x[0] + pt.x[1] + pt.x[2] + pt.x[3])& 1);
	}

	template<class T>
	inline	auto	PointParity(Coord<Colored> pt, Lattice<T,Colored> &myLat) {

		ParityType  pBase = (ParityType)  (( pt.x[0]       +  pt.x[1]      +  pt.x[2]      +  pt.x[3]     ) & 1);
		Parity2Type p2    = (Parity2Type) ((((pt.x[0] & 2) + (pt.x[1] & 2) + (pt.x[2] & 2) + (pt.x[3] & 2)) & 2) >> 1);

		return	(ParityColor) (p2 + ((pt.x[1]&1)<<1) + ((pt.x[2]&1)<<2) + ((pt.x[3]&1)<<3) + (pBase<<4));
	}

	template<class T, ColorKind cOrd>
	inline	auto	PointParity (size_t idx, Lattice<T,cOrd> &myLat) {
		return	PointParity(Coord<cOrd>(idx, myLat.vLength()), myLat);
	}

	template<class T, ColorKind cOrd>
	inline	auto	PointParity (Point pt, Lattice<T,cOrd> &myLat) {
		return	PointParity(pt.idx, myLat);
	}
/*

	template<class T>
	inline	ParityColor	PointColor (Coord pt, Lattice<T> &myLat) {

		ParityType  pBase = (ParityType)  ((  pt.x[0]      +  pt.x[1]      +  pt.x[2]      +  pt.x[3]     ) & 1);
		Parity2Type p2    = (Parity2Type) ((((pt.x[0] & 2) + (pt.x[1] & 2) + (pt.x[2] & 2) + (pt.x[3] & 2)) & 2 ) >> 1);

		return	(ParityColor) (p2 + ((pt.x[1]&1)<<1) + ((pt.x[2]&1)<<2) + ((pt.x[3]&1)<<3) + (pBase<<4));
	}

	template<class T>
	inline	ParityColor	PointColor (size_t idx, Lattice<T> &myLat) {
		return	PointColor(Coord(idx, myLat.vLength()), myLat);
	}

	template<class T>
	inline	ParityColor	PointColor (Point pt, Lattice<T> &myLat) {
		return	PointColor(pt.idx, myLat);
	}

	template<class T>
	inline	ParityType	PointParity (Coord pt, Lattice<T> &myLat) {

		return	(ParityType)  ((pt.x[0] + pt.x[1] + pt.x[2] + pt.x[3])& 1);
	}

	template<class T>
	inline	ParityType	PointParity (size_t idx, Lattice<T> &myLat) {
		return	PointParity(Coord(idx, myLat.vLength()), myLat);
	}

	template<class T>
	inline	ParityType	PointParity (Point pt, Lattice<T> &myLat) {
		return	PointParity(pt.idx, myLat);
	}
*/
#endif
