#ifndef	__SU2_CLASS
	#define	__SU2_CLASS

	#include "random/random.h"

	template<typename Float>
	class	Su2	{
		private:

		Float	a[4];

		public:

		Su2	() 					 { a = { 1., 0., 0., 0. }; }
		Su2	(Float a0, Float a1, Float a2, Float a3) { a[0] = a0;   a[1] = a1;   a[2] = a2;   a[3] = a3;   }
		Su2	(Float *b)				 { a[0] = b[0]; a[1] = b[1]; a[2] = b[2]; a[3] = b[3]; }


		/*	No Vectorization (GPU)	*/

		Su2	operator* (Su2 &b) {
			Su2	tmp;

			tmp.a[0] = a[0]*b.a[0] - a[1]*b.a[1] - a[2]*b.a[2] - a[3]*b.a[3];
			tmp.a[1] = a[0]*b.a[1] + a[1]*b.a[0] - a[2]*b.a[3] + a[3]*b.a[2];
			tmp.a[2] = a[0]*b.a[2] + a[1]*b.a[3] + a[2]*b.a[0] - a[3]*b.a[1];
			tmp.a[3] = a[0]*b.a[3] - a[1]*b.a[2] + a[2]*b.a[1] + a[3]*b.a[0];

			return	tmp;
		}

		Su2&	operator*=(Su2 &b) {
			return	(*this * b);
		}

		Su2	operator+ (Su2 &b) {
			Su2	tmp;

			tmp.a[0] = a[0] + b.a[0];
			tmp.a[1] = a[0] + b.a[1];
			tmp.a[2] = a[0] + b.a[2];
			tmp.a[3] = a[0] + b.a[3];

			return	tmp;
		}

		Su2&	operator+=(Su2 &b) {
			return	(*this + b);
		}

		Su2	operator- (Su2 &b) {
			Su2	tmp;

			tmp.a[0] = a[0] - b.a[0];
			tmp.a[1] = a[0] - b.a[1];
			tmp.a[2] = a[0] - b.a[2];
			tmp.a[3] = a[0] - b.a[3];

			return	tmp;
		}

		Su2&	operator-=(Su2 &b) {
			return	(*this - b);
		}

		Su2	operator* (Float &b) {
			Su2	tmp;

			tmp.a[0] = a[0]*b;
			tmp.a[1] = a[1]*b;
			tmp.a[2] = a[2]*b;
			tmp.a[3] = a[3]*b;

			return	tmp;
		}

		Su2&	operator*=(Float &b) {
			return	(*this * b);
		}

		Float	Norm() {
			return	a[0]*a[0] + a[1]*a[1] + a[2]*a[2] + a[3]*a[3];
		}

		Float	Trace() {
			return	((Float) 2.)*a[0];
		}

		Su2&	Project() {
			return	(*this)/sqrt(Norm());
		}


		Su2&	GenHeat(const Float e2) {

			Float r, s, z;
			const  Float oDet = sqrt(Norm());
			const  Float e    = e2*oDet;
			const  Float b    = exp(((Float) -2.) * e);

			do {
				r    = Su2Rand::genRand();
				s    = (((Float) 1.) - r) * b + r;
				a[0] = ((Float) 1.) + log(s) / e;

				r    = Su2Rand::genRand();
				s    = r*r;
				z    = ((Float) 1.) - a[0]*a[0];
			}	while (s > z);

			s    = sqrt(z);
			r    = Su2Rand::genRand();

			a[3] = (((Float) 2.)*r - ((Float) 1.))*s;

			s    = sqrt(z - a[3]*a[3]);
			r    = ((Float) 2.)*((Float) M_PI)*Su2Rand::genRand();

			a[1] = s*cos(r);
			a[2] = s*sin(r);

			a[0] /= oDet;
			a[1] /= oDet;
			a[2] /= oDet;
			a[3] /= oDet;

			return	(*this);
		}

		Su2&	Pert(const Float eps) {
			#pragma unroll
			for (int i=1; i<4; i++)
				a[i] = ((Float) 2.)*Su2Rand::genRand() - ((Float) 1.);
			a[0] = 1. - sqrt(a[1]*a[1] + a[2]*a[2] + a[3]*a[3]);

			return	(*this);
		}
	};
#endif
