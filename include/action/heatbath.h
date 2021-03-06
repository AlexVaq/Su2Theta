#ifndef	__SU2_HEAT
	#define	__SU2_HEAT

	#include <utils/utils.h>
	#include <action/action.h>
	#include <lattice/lattice.h>

	namespace	Su2Action {

		template<class T, ColorKind cOrd>
		class	HeatBath : public Tunable {
			private:

			Action<T, cOrd> &myAct;

			public:

				HeatBath(Action<T, cOrd> &myAct) : myAct(myAct) {
				InitBlockSize(myAct.Latt().vLength());
				SetName  (std::string("HeatBath"));
				SetColor (std::string(cOrd == Su2Enum::EvenOdd ? "EO" : "P32"));
				SetPrec  (sizeof(typename T::sData) == 4 ? "float" : "double");
				SetVec   (T::sWide != 1 ? Su2Enum::SystemVec : "None");
				SetVolume(myAct.Latt().SLength(), myAct.Latt().TLength());
			}

			void	operator()(const int nSteps = 1)	{
				Su2Prof::Profiler &prof = Su2Prof::getProfiler(ProfHB);
				prof.start();
/*
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
*/
				/*	Test	*/
				size_t cVol = myAct.Latt().oVol()/cOrd;

				for (int i=0; i<nSteps; i++)
				  for (int mu=0; mu<4; mu++)
				    for (int pty=0; pty<cOrd; pty++) {
				      #pragma omp parallel for schedule(static)
				      for (uint idx=0; idx<cVol; idx++) {
					size_t fIdx = cVol*pty + idx;

					myAct.Latt().Insert(std::forward<T>(myAct(fIdx, mu).GenHeat(myAct.Beta())), fIdx, mu);
				      }
				    }
				/*	End test	*/
/*				for (int i=0; i<nSteps; i++)
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
						      auto pPty = PointParity(cPt, myAct.Latt());

						      if (pty != pPty)
							continue;

						      auto idx = cPt.Index(myAct.Latt().vLength());

						      myAct.Latt().Insert(std::forward<T>(myAct(cPt, mu).GenHeat(myAct.Beta())),       idx, mu);
						    }
						  }
						}
					      }
					}
*/
#ifndef	__INTEL_COMPILER
				#undef	Bt
				#undef	Bz
				#undef	By
				#undef	Bx
#endif
				prof.stop();
				double myGFlops = (cOrd == Colored ? 9756.0 : 1500.0) * (myAct.Latt().Volume()* nSteps)*1e-9;
				double myGBytes = (cOrd == Colored ?  436.0 :   76.0) * (myAct.Latt().oVol()  * nSteps)*sizeof(T)/1073741824.0;
				add(myGFlops, myGBytes); 
				prof.add(FullName(), myGFlops, myGBytes);

				LogMsg  (VerbHigh, "%s reporting %lf GFlops %lf GBytes", FullName().c_str(), prof.Prof()[FullName()].GFlops(), prof.Prof()[FullName()].GBytes());
			}

			void	SaveState()	override	{ myAct.Latt().SaveState();    }
			void	RestoreState()	override	{ myAct.Latt().RestoreState(); }
			inline void	Run()	override	{ (*this)(); }
			inline void	Reset()	override	{ Su2Prof::getProfiler(ProfHB).reset(FullName()); }
		};

		template<class T, ColorKind cOrd>
		class	OverRelax : public Tunable {
			private:

			Action<T, cOrd> &myAct;

			public:

				OverRelax(Action<T, cOrd> &myAct) : myAct(myAct) {
				InitBlockSize(myAct.Latt().vLength());
				SetName  (std::string("OverRelax"));
				SetColor (std::string(cOrd == Su2Enum::EvenOdd ? "EO" : "P32"));
				SetPrec  (sizeof(typename T::sData) == 4 ? "float" : "double");
				SetVec   (T::sWide != 1 ? Su2Enum::SystemVec : "None");
				SetVolume(myAct.Latt().SLength(), myAct.Latt().TLength());
			}

			void	operator()(const int nSteps = 1) {
				Su2Prof::Profiler &prof = Su2Prof::getProfiler(ProfOvR);
				prof.start();
/*
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
*/
				/*	Test	*/
				size_t cVol = myAct.Latt().oVol()/cOrd;

				for (int i=0; i<nSteps; i++)
				  for (int mu=0; mu<4; mu++)
				    for (int pty=0; pty<cOrd; pty++) {
				      #pragma omp parallel for schedule(static)
				      for (uint idx=0; idx<cVol; idx++) {
					size_t fIdx = cVol*pty + idx;

					myAct.Latt().Insert(std::forward<T>(myAct.Latt().Data(fIdx, mu).Reflect(std::forward<T>(myAct(fIdx, mu)))), fIdx, mu);
				      }
				    }
				/*	End test	*/
/*				for (int i=0; i<nSteps; i++)
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
						      auto idx = cPt.IdxPt(myAct.Latt().vLength());

						      if (pty != PointParity(cPt, myAct.Latt()))
							continue;

						      myAct.Latt().Insert(std::forward<T>(myAct.Latt().Data(idx, mu).Reflect(std::forward<T>(myAct(cPt, mu)))), idx.idx, mu);
						    }
						  }
						}
					      }
					    }
*/
#ifndef	__INTEL_COMPILER
				#undef	Bt
				#undef	Bz
				#undef	By
				#undef	Bx
#endif
				prof.stop();
				double myGFlops = (cOrd == Colored ? 9720.0 : 1460.0) * (myAct.Latt().Volume()* nSteps)*1e-9;
				double myGBytes = (cOrd == Colored ?  440.0 :   80.0) * (myAct.Latt().oVol()  * nSteps)*sizeof(T)/1073741824.0;
				add(myGFlops, myGBytes); 
				prof.add(FullName(), myGFlops, myGBytes);

				LogMsg  (VerbHigh, "%s reporting %lf GFlops %lf GBytes", FullName().c_str(), prof.Prof()[FullName()].GFlops(), prof.Prof()[FullName()].GBytes());
			}

			void	SaveState()	override	{ myAct.Latt().SaveState();    }
			void	RestoreState()	override	{ myAct.Latt().RestoreState(); }
			inline void	Run()	override	{ (*this)(); }
			inline void	Reset()	override	{ Su2Prof::getProfiler(ProfOvR).reset(FullName()); }
		};
	}
#endif
