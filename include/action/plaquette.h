#ifndef	__SU2_PLAQ
	#define	__SU2_PLAQ

	#include <utils/utils.h>
	#include <lattice/lattice.h>

	using namespace Simd;
	using namespace std;

	namespace	Su2Action {

		template<class T, ColorKind cOrd>
		class	Plaquette : public Tunable {
			private:

			Lattice<T,cOrd> &myLat;

			public:

				Plaquette(Lattice<T,cOrd> &myLat) : myLat(myLat) {
				InitBlockSize(myLat.vLength());
				std::string sName = std::to_string(T::sWide) + std::string(" ") + std::to_string(myLat.SLength()) + std::string("x") + std::to_string(myLat.TLength());
				//SetName(std::string("Plaquette EO\t")  + sName);
				SetName  (std::string("Plaquette"));
				SetColor (std::string(cOrd == Su2Enum::EvenOdd ? "EO" : "P32"));
				SetPrec  (sizeof(typename T::sData) == 4 ? "float" : "double");
				SetVec   (T::sWide != 1 ? Su2Enum::SystemVec : "None");
				SetVolume(myLat.SLength(), myLat.TLength());
			}

			double	operator()()	{
				Su2Prof::Profiler &prof = Su2Prof::getProfiler(ProfPlaq);
				prof.start();

				double	plq = 0.;

				int nBt = myLat.vtLength()/BlockT();
				int nBz = myLat.vzLength()/BlockZ();
				int nBy = myLat.vyLength()/BlockY();
				int nBx = myLat.vxLength()/BlockX();

#ifdef	__INTEL_COMPILER
				const int Bt  = BlockT();
				const int Bz  = BlockZ();
				const int By  = BlockY();
				const int Bx  = BlockX();
#else
				#define	Bt	BlockT()
				#define	Bz	BlockZ()
				#define	By	BlockY()
				#define	Bx	BlockX()
#endif
				for (int bt=0; bt<nBt; bt++)
				  for (int bz=0; bz<nBz; bz++)
				    for (int by=0; by<nBy; by++)
				      for (int bx=0; bx<nBx; bx++) {
					double ret = 0.;
					#pragma omp parallel for collapse(4) reduction(+:ret) schedule(static)//guided)
					for (uint ct=0; ct<Bt; ct++) {
					  for (uint cz=0; cz<Bz; cz++) {
					    for (uint cy=0; cy<By; cy++) {
					      for (uint cx=0; cx<Bx; cx++) {
						auto mPt = Coord<cOrd>(cx + Bx*bx, cy + By*by, cz + Bz*bz, ct + Bt*bt);
						auto mIx = mPt.IdxPt(myLat.vLength());

						for (int mu=0; mu<3; mu++)
						  for (int nu=mu+1; nu<4; nu++) {
						    auto nPt = myLat.nextPoint(mPt, mu);
						    auto oPt = myLat.nextPoint(mPt, nu);

						    auto tmp = myLat.Data(mIx, mu);
						    tmp *=  myLat.Data(nPt, nu);
						    tmp *= !myLat.Data(oPt, mu);
						    tmp *= !myLat.Data(mIx, nu);

						    ret += tmp.SuperTrace();
						}
					      }
					    }
					  }
				        }

					plq += ret;
				      }
#ifndef	__INTEL_COMPILER
				#undef	Bt
				#undef	Bz
				#undef	By
				#undef	Bx
#endif
				prof.stop();
				double myGFlops = 517.0 * myLat.Volume() * 1e-9;
				double myGBytes =  24.0 * myLat.oVol()   * sizeof(T) / 1073741824.0;
				add(myGFlops, myGBytes); 
				prof.add(FullName(), myGFlops, myGBytes);

				LogMsg  (VerbHigh, "%s reporting %lf GFlops %lf GBytes", FullName().c_str(), prof.Prof()[FullName()].GFlops(), prof.Prof()[FullName()].GBytes());
				return	plq*cfPlq/(6.*myLat.Volume());
			}

			void	SaveState()	override	{}
			void	RestoreState()	override	{}
			inline	void	Run()	override	{ (*this)(); }
			inline void	Reset() override	{ Su2Prof::getProfiler(ProfPlaq).reset(FullName()); }
		};
	}
#endif
