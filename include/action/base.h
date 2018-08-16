#ifndef	__SU2_BASE_ACTION
	#define	__SU2_BASE_ACTION

	#include <enumFields.h>
	#include <lattice/lattice.h>

	namespace Su2Action {

		constexpr double cfWilson   = 0.5/6;
		constexpr double cfImproved = 0.;

		template<class T>
		class	Action {
			private:

			double	beta;
			double	theta;

			const bool impAction;

			Lattice<T> &lat;

			public:

				Action (Lattice<T> &lat, double beta, double theta, bool imp) : lat(lat), beta(beta), theta(theta), impAction(imp) {};

			double		Beta () { return beta;  }
			double		Theta() { return theta; }
			Lattice<T>&	Latt()  { return lat; }

			// Con beta
			double	operator()() {
				return	allPts()*beta;
			}

			inline double	operator()(size_t mPt) {
				return	(*this)[mPt]*beta;
			}

			// Sin beta
			inline double	allPts() {
				double	total = 0.;

				// Divide into blocks
				#pragma omp parallel for reduction(+:total) schedule(static)
				for (size_t i=0; i<lat.Volume(); i++)
					total += (*this)[i];

				return	total;
			}

			inline T	operator()(size_t mPt, int mu) {

				T	staple(0.), iSp(0.);

				for (int dir=mu+1; dir<mu+4; dir++) {
					int nu = dir%4;

					auto nPt = lat.nextPoint(mPt, mu);
					auto oPt = lat.nextPoint(mPt, nu);

					auto tmp = lat.Data(nPt, nu);
					tmp *= !lat.Data(oPt, mu);
					tmp *= !lat.Data(mPt, nu);

					staple += tmp;

					auto pPt = lat.nextPoint(mPt, (nu+4));
					auto qPt = lat.nextPoint(nPt, (nu+4));

					tmp  = !lat.Data(qPt, nu);
					tmp *= !lat.Data(pPt, mu);
					tmp *=  lat.Data(pPt, nu);

					staple += tmp;
				}

				if (impAction) {
					for (int mu=0; mu<3; mu++) {
						for (int nu=mu+1; nu<4; nu++) {
							auto nPt = lat.nextPoint(mPt, mu);
							auto oPt = lat.nextPoint(nPt, mu);
							auto qPt = lat.nextPoint(mPt, nu);
							auto pPt = lat.nextPoint(nPt, nu);

							auto tmp = lat.Data(mPt, mu);
							tmp *= lat.Data(nPt, mu);
							tmp *= lat.Data(oPt, nu);
							tmp *= !lat.Data(pPt, mu);
							tmp *= !lat.Data(qPt, mu);
							tmp *= !lat.Data(mPt, nu);

							iSp += tmp;

							oPt  = lat.nextPoint(nPt, nu);
							pPt  = lat.nextPoint(qPt, nu);

							tmp  = lat.Data(mPt, mu);
							tmp *= lat.Data(nPt, nu);
							tmp *= lat.Data(oPt, nu);
							tmp *= !lat.Data(pPt, mu);
							tmp *= !lat.Data(qPt, nu);
							tmp *= !lat.Data(mPt, nu);

							iSp += tmp;
						}
					}

					staple += iSp*(cfImproved/cfWilson);
				}

				return	staple;
			}

			inline double	operator[](size_t mPt) {

				double	ret = 0.;
				T	tmp;

				for (int mu=0; mu<3; mu++) {
					for (int nu=mu+1; nu<4; nu++) {
						auto nPt = lat.nextPoint(mPt, mu);
						auto oPt = lat.nextPoint(mPt, nu);

						tmp  =  lat.Data(mPt, mu);
						tmp *=  lat.Data(nPt, nu);
						tmp *= !lat.Data(oPt, mu);
						tmp *= !lat.Data(mPt, nu);

						ret += tmp.SuperTrace();
					}
				}

				ret *= cfWilson;

				if (impAction) {
					double	iRt = 0.;

					for (int mu=0; mu<3; mu++) {
						for (int nu=mu+1; nu<4; nu++) {
							auto nPt = lat.nextPoint(mPt, mu);
							auto oPt = lat.nextPoint(nPt, mu);
							auto qPt = lat.nextPoint(mPt, nu);
							auto pPt = lat.nextPoint(nPt, nu);

							tmp  =  lat.Data(mPt, mu);
							tmp *=  lat.Data(nPt, mu);
							tmp *=  lat.Data(oPt, nu);
							tmp *= !lat.Data(pPt, mu);
							tmp *= !lat.Data(qPt, mu);
							tmp *= !lat.Data(mPt, nu);

							iRt += tmp.SuperTrace();

							oPt  = lat.nextPoint(nPt, nu);
							pPt  = lat.nextPoint(qPt, nu);

							tmp  =  lat.Data(mPt, mu);
							tmp *=  lat.Data(nPt, nu);
							tmp *=  lat.Data(oPt, nu);
							tmp *= !lat.Data(pPt, mu);
							tmp *= !lat.Data(qPt, nu);
							tmp *= !lat.Data(mPt, nu);

							iRt += tmp.SuperTrace();
						}
					}

					ret += iRt*cfImproved;
				}

				return	ret;
			}
		};
	}
#endif

