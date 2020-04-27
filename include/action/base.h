#ifndef	__SU2_BASE_ACTION
	#define	__SU2_BASE_ACTION

	#include <enumFields.h>
	#include <utils/utils.h>
	#include <lattice/lattice.h>

	namespace Su2Action {

		constexpr double cfPlq     =  0.5;
		constexpr double cfImpPlq  =  5./3.;
		constexpr double cfImpRect = -1./12.;

		template<class T, ColorKind cOrd>
		class	Action  : public Tunable {
			private:

			double	beta;
			double	theta;

			Lattice<T, cOrd> &lat;

			public:

			typedef T Data;

				Action (Lattice<T, cOrd> &lat, double beta, double theta) : lat(lat), beta(beta), theta(theta) {
				InitBlockSize(lat.vLength());
				SetName  (std::string("Action"));
				SetColor (std::string(cOrd == Su2Enum::EvenOdd ? "EO" : "P32"));
				SetPrec  (sizeof(typename T::sData) == 4 ? "float" : "double");
				SetVec   (T::sWide != 1 ? Su2Enum::SystemVec : "None");
				SetVolume(lat.SLength(), lat.TLength());
				/*
				std::string sName = std::to_string(T::sWide) + std::string(" ") + std::to_string(lat.SLength()) + std::string("x") + std::to_string(lat.TLength());
				if (cOrd == Su2Enum::EvenOdd)
					SetName(std::string("Action EO\t")  + sName);
				else
					SetName(std::string("Action P32\t") + sName);
				*/
			};

			inline	double		Beta () { return beta;  }
			inline	double		Theta() { return theta; }
			inline	Lattice<T,cOrd>&Latt()  { return lat; }

			inline	constexpr ColorKind	ParityLevel() const { return cOrd; }

			// Con beta
//			double	operator()() {
//				return	allPts()*beta;
//			}

//			inline	double	operator()(size_t mPt) {
//				return	(*this)[mPt]*beta;
//			}

			// Sin beta
			inline	double	allPts() {

				Su2Prof::Profiler &prof = Su2Prof::getProfiler(ProfAction);
				prof.start();

				double	total = 0.;

				// Divide into blocks
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
				const int nBt = Latt().vtLength()/BlockT();
				const int nBz = Latt().vzLength()/BlockZ();
				const int nBy = Latt().vyLength()/BlockY();
				const int nBx = Latt().vxLength()/BlockX();

				for (int bt=0; bt<nBt; bt++)
				  for (int bz=0; bz<nBz; bz++)
				    for (int by=0; by<nBy; by++)
				      for (int bx=0; bx<nBx; bx++) {
					double partial = 0.;

					#pragma omp parallel for reduction(+:partial) collapse(4) schedule(static)
					for (uint ct=0; ct<Bt; ct++) {
					  for (uint cz=0; cz<Bz; cz++) {
					    for (uint cy=0; cy<By; cy++) {
					      for (uint cx=0; cx<Bx; cx++) {
						auto iPt = Coord<cOrd>(cx + BlockX()*bx, cy + BlockY()*by, cz + BlockZ()*bz, ct + BlockT()*bt);
						partial += (*this)(iPt);
					      }
					    }
					  }
					}
					total += partial;
				      }
#ifndef	__INTEL_COMPILER
				#undef	Bt
				#undef	Bz
				#undef	By
				#undef	Bx
#endif
				prof.stop();
				double myGFlops = (cOrd == Colored ? 2197.0 : 517.0) * Latt().Volume() * 1e-9;
				double myGBytes = (cOrd == Colored ?   96.0 :  24.0) * Latt().oVol()   * sizeof(T) / 1073741824.0;

				add(myGFlops, myGBytes); 
				prof.add(FullName(), myGFlops, myGBytes);

				LogMsg  (VerbHigh, "%s reporting %lf GFlops %lf GBytes", FullName().c_str(), prof.Prof()[FullName()].GFlops(), prof.Prof()[FullName()].GBytes());

				return	total;
			}

			inline	T	operator()(size_t mPt, int mu) {
				return (*this)(Coord<cOrd>(mPt, lat.vLength()), mu);
			}

			inline	T	operator()(Coord<cOrd> mPt, int mu) {

				T	staple(0.), iSp(0.);

				for (int dir=mu+1; dir<mu+4; dir++) {
					Coord<cOrd> pPt, qPt;
					int nu = dir%4;

					auto nPt = lat.nextPoint(mPt, mu);
					auto oPt = lat.nextPoint(mPt, nu);
					auto mIx = mPt.IdxPt(lat.vLength());
					auto nIx = nPt.IdxPt(lat.vLength());
					auto oIx = oPt.IdxPt(lat.vLength());
					auto pIx = oIx;
					auto qIx = oIx;

					auto tmp = lat.Data(nIx, nu);
					tmp *= !lat.Data(oIx, mu);
					tmp *= !lat.Data(mIx, nu);

					staple += tmp;

					if (cOrd == Su2Enum::Colored) {
						qPt = lat.nextPoint(nPt, mu);
						pPt = lat.nextPoint(nPt, nu);
						pIx = pPt.IdxPt(lat.vLength());

						tmp  =  lat.Data(nIx, mu);
						tmp *=  lat.Data(qPt, nu);
						tmp *= !lat.Data(pIx, mu);
						tmp *= !lat.Data(oIx, mu);
						tmp *= !lat.Data(mIx, nu);

						iSp += tmp;

						qPt  =  lat.nextPoint(oPt, nu);

						tmp  =  lat.Data(nIx, nu);
						tmp *=  lat.Data(pIx, nu);
						tmp *= !lat.Data(qPt, mu);
						tmp *= !lat.Data(oIx, nu);
						tmp *= !lat.Data(mIx, nu);

						iSp += tmp;

						qPt = lat.nextPoint(oPt, mu+4);
						pPt = lat.nextPoint(mPt, mu+4);
						qIx = qPt.IdxPt(lat.vLength());
						pIx = pPt.IdxPt(lat.vLength());

						tmp  =  lat.Data(nIx, nu);
						tmp *= !lat.Data(oIx, mu);
						tmp *= !lat.Data(qPt, mu);
						tmp *= !lat.Data(pIx, nu);
						tmp *=  lat.Data(pIx, mu);

						iSp += tmp;
					}

					// Opposite direction
					auto rPt = lat.nextPoint(mPt, (nu+4));
					oPt = lat.nextPoint(nPt, (nu+4));
					auto rIx = rPt.IdxPt(lat.vLength());
					oIx = oPt.IdxPt(lat.vLength());

					tmp  = !lat.Data(oIx, nu);
					tmp *= !lat.Data(rIx, mu);
					tmp *=  lat.Data(rIx, nu);

					staple += tmp;

					if (cOrd == Su2Enum::Colored) {
						qPt = lat.nextPoint(rPt, (mu+4));
						qIx = qPt.IdxPt(lat.vLength());

						tmp  = !lat.Data(oIx, nu);
						tmp *= !lat.Data(rIx, mu);
						tmp *= !lat.Data(qIx, mu);
						tmp *=  lat.Data(qIx, nu);
						tmp *=  lat.Data(pIx, mu);

						iSp += tmp;

						pPt = lat.nextPoint(oPt, nu+4);
						qPt = lat.nextPoint(rPt, nu+4);
						qIx = qPt.IdxPt(lat.vLength());

						tmp  = !lat.Data(oIx, nu);
						tmp *= !lat.Data(pPt, nu);
						tmp *= !lat.Data(qIx, mu);
						tmp *=  lat.Data(qIx, nu);
						tmp *=  lat.Data(rIx, nu);

						iSp += tmp;

						qPt = lat.nextPoint(oPt, mu);

						tmp  =  lat.Data(nIx, mu);
						tmp *= !lat.Data(qPt, nu);
						tmp *= !lat.Data(oIx, mu);
						tmp *= !lat.Data(rIx, mu);
						tmp *=  lat.Data(rIx, nu);

						iSp += tmp;
					}
				}

				return	(cOrd == Colored ? (staple*cfImpPlq + iSp*cfImpRect) : staple);
			}

			inline double	operator()(Coord<cOrd> mPt) {

				double	ret = 0., iRt = 0.;
				T	tmp;

				for (int mu=0; mu<3; mu++) {
					for (int nu=mu+1; nu<4; nu++) {
						auto nPt = lat.nextPoint(mPt, mu);
						auto oPt = lat.nextPoint(mPt, nu);

						auto mIx = mPt.IdxPt(lat.vLength());
						auto nIx = nPt.IdxPt(lat.vLength());
						auto oIx = oPt.IdxPt(lat.vLength());

						tmp  =  lat.Data(mIx, mu);
						tmp *=  lat.Data(nIx, nu);
						tmp *= !lat.Data(oIx, mu);
						tmp *= !lat.Data(mIx, nu);

						ret += tmp.SuperTrace();

						if (cOrd == Su2Enum::Colored) {
							auto pPt = lat.nextPoint(nPt, mu);
							auto qPt = lat.nextPoint(nPt, nu);
							auto rPt = lat.nextPoint(oPt, nu);

							auto qIx = qPt.IdxPt(lat.vLength());

							tmp  =  lat.Data(mIx, mu);
							tmp *=  lat.Data(nIx, mu);
							tmp *=  lat.Data(pPt, nu);
							tmp *= !lat.Data(qIx, mu);
							tmp *= !lat.Data(oIx, mu);
							tmp *= !lat.Data(mIx, nu);

							iRt += tmp.SuperTrace();

							tmp  =  lat.Data(mIx, mu);
							tmp *=  lat.Data(nIx, nu);
							tmp *=  lat.Data(qIx, nu);
							tmp *= !lat.Data(rPt, mu);
							tmp *= !lat.Data(oIx, nu);
							tmp *= !lat.Data(mIx, nu);

							iRt += tmp.SuperTrace();
						}
					}
				}

				return	(cOrd == Colored ? (ret*(cfImpPlq*cfPlq) + iRt*(cfImpRect*cfPlq)) : ret*cfPlq);
			}

			inline double	operator()(size_t mPt) {
				return	(*this)(Coord<cOrd>(mPt, lat.vLength()));

			}

			void	SaveState()	override	{}
			void	RestoreState()	override	{}
			inline	void	Run()	override	{ (*this).allPts(); }
			inline void	Reset()	override	{ Su2Prof::getProfiler(ProfAction).reset(FullName()); }
		};

	}
#endif

