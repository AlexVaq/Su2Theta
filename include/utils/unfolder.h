#ifndef	__UTIL_UNFOLD
	#define	__UTIL_UNFOLD

	#include <utils/utils.h>
	#include <lattice/lattice.h>

	template<typename vT, ColorKind cOrd>
	class	Unfolder : public Tunable {
		private:

		Lattice<vT,cOrd>	&vLat;

		public:

		typedef typename vT::sClass T;

			Unfolder (Lattice<vT,cOrd> &lat) : vLat(lat) {
			std::string sName = std::to_string(T::sWide) + std::string(" ") + std::to_string(lat.SLength()) + std::string("x") + std::to_string(lat.TLength());
			SetName(std::string("Folder EO\t")  + sName);
		}

		Lattice<T,cOrd>	operator()() {
			Su2Prof::Profiler &prof = Su2Prof::getProfiler(ProfFold);
			prof.start();

			auto	vLx = vLat.vxLength();
			auto	vLy = vLat.vyLength();
			auto	vLz = vLat.vzLength();
			auto	vLt = vLat.vtLength();
			auto	vS  = vLx*vLy;
			auto	vV  = vS *vLz;

			auto	Ls  = vLat.SLength();
			auto	Lt  = vLat.SLength();
			auto	S   = Ls*Ls;
			auto	Vs  = vLat.SVol();

			Lattice<T,cOrd> sLat(vLat.SLength(), vLat.TLength());

			for (size_t mu=0; mu<4; mu++) {
				#pragma omp parallel for schedule(static) 
				for (size_t i=0; i<vLat.oVol(); i++) {
				  auto U  = vLat.Data(i,mu);

				  for (int wx=0; wx<vT::xWide; wx++)
				    for (int wy=0; wy<vT::yWide; wy++)
				      for (int wz=0; wz<vT::zWide; wz++)
					for (int wt=0; wt<vT::tWide; wt++) {
					  int lane = (wt + vT::tWide*(wz + vT::zWide*(wy + vT::yWide*wx)));
					  Coord<cOrd> vPt(i, vLat.vLength());

					  vPt.x[0] += wx*vLx;	vPt.x[1] += wy*vLy;
					  vPt.x[2] += wz*vLz;	vPt.x[3] += wt*vLt;

					  sLat.Insert(U.Extract(lane), vPt.Index(sLat.vLength()), mu);
					}
				}
			}

			prof.stop();
			double myGBytes = 8.0 * vLat.oVol() * sizeof(T) / 1073741824.0;
			add(0., myGBytes); 
			prof.add(Name(), 0., myGBytes);

			LogMsg  (VerbHigh, "%s reporting %lf GFlops %lf GBytes", Name().c_str(), prof.Prof()[Name()].GFlops(), prof.Prof()[Name()].GBytes());

			return	sLat;
		}

		void	SaveState()	override	{}
		void	RestoreState()	override	{}
		inline void	Run()	override	{ (*this)(); }
		inline void	Reset() override	{ Su2Prof::getProfiler(ProfFold).reset(Name()); }
	};
#endif
