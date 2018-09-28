#ifndef	__SU2_METRO
	#define	__SU2_METRO

	#include <action/action.h>
	#include <random/random.h>
	#include <lattice/lattice.h>

	//constexpr int nTries = 4;

	using namespace Simd;
	using namespace std;

	namespace	Su2Action {

		template<class T, ColorKind cOrd>
		class	Metropolis : public Tunable {
			private:

			Action<T,cOrd> &myAct;

			public:

				Metropolis(Action<T,cOrd> &myAct) : myAct(myAct) {
					InitBlockSize(myAct.Latt().vLength());
					std::string sName = std::to_string(T::sWide) + std::string(" ") + std::to_string(myAct.Latt().SLength()) + std::string("x") + std::to_string(myAct.Latt().TLength());
					if (cOrd == Su2Enum::EvenOdd)
						SetName(std::string("Metro EO \t")  + sName);
					else
						SetName(std::string("Metro P32\t") + sName);
			}

			void	operator()(const int nTries = 4)	{
				Su2Prof::Profiler &prof = Su2Prof::getProfiler(ProfMetro);
				prof.start();

				int nBt = myAct.Latt().vtLength()/BlockT();
				int nBz = myAct.Latt().vzLength()/BlockZ();
				int nBy = myAct.Latt().vyLength()/BlockY();
				int nBx = myAct.Latt().vxLength()/BlockX();

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
				      for (int bx=0; bx<nBx; bx++)
					for (int mu=0; mu<4; mu++)
					  for (int pty=0; pty<cOrd; pty++) {
					    #pragma omp parallel for collapse(4) schedule(static)//guided)
					    for (uint ct=0; ct<Bt; ct++) {
					      for (uint cz=0; cz<Bz; cz++) {
						for (uint cy=0; cy<By; cy++) {
						  for (uint cx=0; cx<Bx; cx++) {
						    auto cPt = Coord<cOrd>(cx + Bx*bx, cy + By*by, cz + Bz*bz, ct + Bt*bt);

						    if (pty != PointParity(cPt, myAct.Latt()))
							continue;

						    auto stpl = myAct(cPt, mu);

						    #pragma unroll
						    for (int j=0; j<nTries; j++) {
							auto link = myAct.Latt().Data(cPt, mu);
							auto nLnk = link.Pert(1.0);
							auto dLnk = ((nLnk - link)*stpl).Trace()*(myAct.Beta()*(cOrd == EvenOdd ? 0.5 : 5./6.));
							//myAct.Latt().InsertMask((Su2Rand::genVRand<typename T::data>() < exp(dLnk)), std::forward<T>(nLnk), cPt.Index(myAct.Latt().vLength()), mu);
							myAct.Latt().InsertMask((Su2Rand::vGenRand<typename T::data>() < exp(dLnk)), std::forward<T>(nLnk), cPt.Index(myAct.Latt().vLength()), mu);
						    }
						  }
						}
					      }
					    }
					}
#ifndef	__INTEL_COMPILER
				#undef	Bt
				#undef	Bz
				#undef	By
				#undef	Bx
#endif
				prof.stop();
				double myGFlops = ((cOrd == Colored ? 9728.0 : 1428.0) + 34.0 * nTries) * myAct.Latt().Volume() * 1e-9;
				double myGBytes = ((cOrd == Colored ?  432.0 :   72.0) +  8.0 * nTries) * myAct.Latt().oVol() * sizeof(T) / 1073741824.0;
				add(myGFlops, myGBytes); 
				prof.add(Name(), myGFlops, myGBytes);

				LogMsg  (VerbHigh, "%s reporting %lf GFlops %lf GBytes", Name().c_str(), prof.Prof()[Name()].GFlops(), prof.Prof()[Name()].GBytes());
			}

			void	SaveState()	override	{ myAct.Latt().SaveState();    }
			void	RestoreState()	override	{ myAct.Latt().RestoreState(); }
			inline	void	Run()	override	{ (*this)(); }
			inline void	Reset() override	{ Su2Prof::getProfiler(ProfMetro).reset(Name()); }
		};
	}
#endif
