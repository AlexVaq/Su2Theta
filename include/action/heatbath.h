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
				/*
				std::string sName = std::to_string(T::sWide) + std::string(" ") + std::to_string(myAct.Latt().SLength()) + std::string("x") + std::to_string(myAct.Latt().TLength());
				if (cOrd == Su2Enum::EvenOdd)
					SetName(std::string("HeatBath EO\t")  + sName);
				else
					SetName(std::string("HeatBath P32\t") + sName);
				*/
			}

			void	operator()(const int nSteps = 1)	{
				Su2Prof::Profiler &prof = Su2Prof::getProfiler(ProfHB);
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
				for (int i=0; i<nSteps; i++)
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

//						      switch (cOrd) {
//							case EvenOdd:
							myAct.Latt().Insert(std::forward<T>(myAct(cPt, mu).GenHeat(myAct.Beta())),       idx, mu);
//							break;

//							case Colored:
//							myAct.Latt().Insert(std::forward<T>(myAct(cPt, mu).GenHeat(myAct.Beta()*5./3.)), idx, mu);
//							break;
  //						      }
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
				/*
				std::string sName = std::to_string(T::sWide) + std::string(" ") + std::to_string(myAct.Latt().SLength()) + std::string("x") + std::to_string(myAct.Latt().TLength());
				if (cOrd == Su2Enum::EvenOdd)
					SetName(std::string("OverRelax EO\t")  + sName);
				else
					SetName(std::string("OverRelax P32\t") + sName);
				*/
			}

			void	operator()(const int nSteps = 1) {
				Su2Prof::Profiler &prof = Su2Prof::getProfiler(ProfOvR);
				prof.start();

				int nBt = 1;//myAct.Latt().vtLength()/BlockT();
				int nBz = 1;//myAct.Latt().vzLength()/BlockZ();
				int nBy = 1;//myAct.Latt().vyLength()/BlockY();
				int nBx = 1;//myAct.Latt().vxLength()/BlockX();

#ifdef	__INTEL_COMPILER
				const int Bt  = 1;//BlockT();
				const int Bz  = 1;//BlockZ();
				const int By  = 1;//BlockY();
				const int Bx  = 1;//BlockX();
#else
				#define	Bt	1 //BlockT()
				#define	Bz	1 //BlockZ()
				#define	By	1 //BlockY()
				#define	Bx	1 //BlockX()
#endif
				for (int i=0; i<nSteps; i++)
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
						    //  auto staple = myAct(cPt, mu);
						    //  auto link   = myAct.Latt().Data(cPt, mu);
						    //  auto nLink  = link.Reflect(std::forward<T>(staple));
						    //  auto diff1  = (staple*(link - nLink)).SuperTrace()*0.5;
						    //  auto diff2  = myAct.allPts();
						    //  diff1 = sqrt(diff1*diff1);
						    //  if (diff1 > 1e-12)
						    //    printf("CACA %le\n", diff1);
						    //  //else
						    //    //printf("%d (%d %d %d %d) %d\n", mu, cPt.x[0], cPt.x[1], cPt.x[2], cPt.x[3], pty);
						    //  myAct.Latt().Insert(std::forward<T>(nLink), cPt.Index(myAct.Latt().vLength()), mu);
						    //  diff2 -= myAct.allPts();
						    //  diff2 = sqrt(diff2*diff2);
						    //  if (diff2 > 1e-12) {
						    //    printf("COCO %le vs %le %d/(%d %d %d %d) - %d\n", diff2, diff1, mu, cPt.x[0], cPt.x[1], cPt.x[2], cPt.x[3], cPt.Index(myAct.Latt().vLength()));
						    //    //auto rLink = nLink - myAct.Latt().Data(cPt, mu);
						    //    //rLink.Print();
						    //  }
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
