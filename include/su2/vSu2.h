#ifndef	__VSU2_CLASS
	#define	__VSU2_CLASS

	#include "random/random.h"
	#include <algorithm>

	template<typename vFloat>
	class	vSu2	{
		private:

		vFloat	a[4];

		public:

		typedef		 vFloat        data;
		typedef	typename vFloat::sData sData;

		static constexpr size_t sWide = sizeof(vFloat)/sizeof(sData);

		vSu2	() 					     { a[0] = vFloat(1.); a[1] = a[2] = a[3] = vFloat(0.); }
		vSu2	(vFloat a0, vFloat a1, vFloat a2, vFloat a3) { a[0] = a0;   a[1] = a1;   a[2] = a2;   a[3] = a3;   }
		vSu2	(vFloat *b)				     { for(int i=0;i<4;i++) a[i] = b[i]; }


		/*	Vectorized datatypes(CPU)	*/

		vSu2	operator* (const vSu2 &b) {
			vSu2	tmp;

			tmp.a[0] = a[0]*b.a[0] - a[1]*b.a[1] - a[2]*b.a[2] - a[3]*b.a[3];
			tmp.a[1] = a[0]*b.a[1] + a[1]*b.a[0] - a[2]*b.a[3] + a[3]*b.a[2];
			tmp.a[2] = a[0]*b.a[2] + a[1]*b.a[3] + a[2]*b.a[0] - a[3]*b.a[1];
			tmp.a[3] = a[0]*b.a[3] - a[1]*b.a[2] + a[2]*b.a[1] + a[3]*b.a[0];

			return	tmp;
		}

		vSu2	operator*=(const vSu2 &b) {
			return	(*this * b);
		}

		vSu2	operator+ (const vSu2 &b) {
			vSu2	tmp;

			tmp.a[0] = a[0] + b.a[0];
			tmp.a[1] = a[0] + b.a[1];
			tmp.a[2] = a[0] + b.a[2];
			tmp.a[3] = a[0] + b.a[3];

			return	tmp;
		}

		vSu2	operator+=(const vSu2 &b) {
			return	(*this + b);
		}

		vSu2	operator- (const vSu2 &b) {
			vSu2	tmp;

			tmp.a[0] = a[0] - b.a[0];
			tmp.a[1] = a[0] - b.a[1];
			tmp.a[2] = a[0] - b.a[2];
			tmp.a[3] = a[0] - b.a[3];

			return	tmp;
		}

		vSu2	operator! () {
			vSu2	tmp;

			tmp.a[1] = -a[1];
			tmp.a[2] = -a[2];
			tmp.a[3] = -a[3];

			return	tmp;
		}

		vSu2	operator-=(const vSu2 &b) {
			return	(*this - b);
		}

		vSu2	operator* (const vFloat &b) {
			vSu2	tmp;

			tmp.a[0] = a[0]*b;
			tmp.a[1] = a[1]*b;
			tmp.a[2] = a[2]*b;
			tmp.a[3] = a[3]*b;

			return	tmp;
		}

		vSu2	operator*=(const vFloat &b) {
			return	(*this * b);
		}

		vFloat	Norm() {
			return	a[0]*a[0] + a[1]*a[1] + a[2]*a[2] + a[3]*a[3];
		}

		vFloat	Trace() {
			return	vFloat(2.)*a[0];
		}

		sData	SuperTrace() {	// Only different for vectorized datatypes
			return	Trace().Sum();
		}

		vSu2	Project() {
			return	(*this)/sqrt(Norm());
		}


		vSu2&	GenHeat(const sData e2) {

			vFloat r, s, z;
			const  vFloat oDet = sqrt(Norm());
			const  vFloat e    = oDet*vFloat(e2);
			const  vFloat b    = exp(e*vFloat(-2.));

			do {
				r    = Su2Rand::genVRand<vFloat>();
				s    = (vFloat(1.) - r) * b + r;
				a[0] = vFloat(1.) + log(s) / e;

				r    = Su2Rand::genVRand<vFloat>();
				s    = r*r;
				z    = vFloat(1.) - a[0]*a[0];
			}	while (s > z);

			s    = sqrt(z);
			r    = Su2Rand::genVRand<vFloat>();

			a[3] = (vFloat(2.)*r - vFloat(1.))*s;

			s    = sqrt(z - a[3]*a[3]);
			r    = vFloat(2.)*vFloat(M_PI)*Su2Rand::genVRand<vFloat>();

			a[1] = s*cos(r);
			a[2] = s*sin(r);

			a[0] /= oDet;
			a[1] /= oDet;
			a[2] /= oDet;
			a[3] /= oDet;

			return	(*this);
		}

		vSu2&	SetRandom(const sData eps = 1.) {
			#pragma unroll
			for (int i=0; i<4; i++) {
				vFloat r = Su2Rand::genVRand<vFloat>();
				a[i] = vFloat(2.*eps)*r - vFloat(eps);
			}

			auto mod = sqrt(vFloat(1.)/(a[0]*a[0] + a[1]*a[1] + a[2]*a[2] + a[3]*a[3]));

			a[0] *= mod;
			a[1] *= mod;
			a[2] *= mod;
			a[3] *= mod;

			return	(*this);
		}

		vSu2&	Pert	(const typename vFloat::sData eps) {
			vSu2<vFloat> tmp;

			tmp.SetRandom(eps);

			(*this) *= tmp;
		}

		vSu2	xPermute() {
			vSu2<vFloat> tmp;
			#pragma unroll
			for (int i=0; i<4; i++)
				tmp.a[i] = a[i].xPermute();
			return	tmp;
		}

		vSu2	yPermute() {
			vSu2<vFloat> tmp;
			#pragma unroll
			for (int i=0; i<4; i++)
				tmp.a[i] = a[i].yPermute();
			return	tmp;
		}

		vSu2	zPermute() {
			vSu2<vFloat> tmp;
			#pragma unroll
			for (int i=0; i<4; i++)
				tmp.a[i] = a[i].zPermute();
			return	tmp;
		}

		vSu2	tPermute() {
			vSu2<vFloat> tmp;
			#pragma unroll
			for (int i=0; i<4; i++)
				tmp.a[i] = a[i].tPermute();
			return	tmp;
		}
	};
#endif
