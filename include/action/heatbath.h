#ifndef	__SU2_HEAT
	#define	__SU2_HEAT

	#include <action/action.h>
	#include <lattice/lattice.h>

	using namespace Su2Action;

	namespace	Su2HB {

		template<class T>
		void	HeatBath (Action<T> &myAct) {
			// Parity loop, for non-improved action
			for (int pty=0; pty<2; pty++) {
				for (int mu=0; mu<4; mu++) {
					#pragma omp parallel for schedule(static)
					for (size_t i=0; i<myAct.Latt().oVol(); i++) {
						auto pPty = PointParity(i, myAct.Latt());

						if (pty != ((ParityType) pPty))
							continue;
						//auto link = myAct(i, mu).GenHeat(myAct.Beta());
						//myAct.Latt().Insert(link, i, mu);
						myAct.Latt().Insert(myAct(i, mu).GenHeat(myAct.Beta()), i, mu);
					}
				}
			}
		}
	}
#endif
