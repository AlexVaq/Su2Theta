#ifndef	__SU2_CLASS
	#define	__SU2_CLASS

	#include "random/random.h"
	#include <algorithm>

	template<typename Float>
	class	Su2	{
		private:

		Float	a[4];

		public:

		typedef	Float  data;
		typedef	Float sData;

		static constexpr size_t sWide = 1;
		static constexpr size_t xWide = 1;
		static constexpr size_t yWide = 1;
		static constexpr size_t zWide = 1;
		static constexpr size_t tWide = 1;

		Su2	() 					 { a[0] = 1.; a[1] = a[2] = a[3] = 0.; }
		Su2	(Float a0, Float a1, Float a2, Float a3) { a[0] = a0;   a[1] = a1;   a[2] = a2;   a[3] = a3;   }
		Su2	(Float *b)				 { std::copy(b, b+4, a); }
		Su2	(Float  c)				 { a[0] = a[1] = a[2] = a[3] = c; }


		/*	No Vectorization (GPU)	*/

		Su2	operator* (const Su2 &b) {
			Su2	tmp;

			tmp.a[0] = a[0]*b.a[0] - a[1]*b.a[1] - a[2]*b.a[2] - a[3]*b.a[3];
			tmp.a[1] = a[0]*b.a[1] + a[1]*b.a[0] - a[2]*b.a[3] + a[3]*b.a[2];
			tmp.a[2] = a[0]*b.a[2] + a[1]*b.a[3] + a[2]*b.a[0] - a[3]*b.a[1];
			tmp.a[3] = a[0]*b.a[3] - a[1]*b.a[2] + a[2]*b.a[1] + a[3]*b.a[0];

			return	tmp;
		}

		Su2	operator*=(const Su2 &b) {
			(*this) = (*this)*b;
			return	(*this);
		}

		Su2	operator+ (const Su2 &b) {
			Su2	tmp;

			tmp.a[0] = a[0] + b.a[0];
			tmp.a[1] = a[0] + b.a[1];
			tmp.a[2] = a[0] + b.a[2];
			tmp.a[3] = a[0] + b.a[3];

			return	tmp;
		}

		Su2	operator+=(const Su2 &b) {
			(*this) = (*this)+b;
			return	(*this);
		}

		Su2	operator- (const Su2 &b) {
			Su2	tmp;

			tmp.a[0] = a[0] - b.a[0];
			tmp.a[1] = a[0] - b.a[1];
			tmp.a[2] = a[0] - b.a[2];
			tmp.a[3] = a[0] - b.a[3];

			return	tmp;
		}

		Su2	operator-=(const Su2 &b) {
			(*this) = (*this)-b;
			return	(*this);
		}

		Su2	operator! () {
			Su2	tmp;

			tmp.a[0] =  a[0];
			tmp.a[1] = -a[1];
			tmp.a[2] = -a[2];
			tmp.a[3] = -a[3];

			return	tmp;
		}

		Su2	operator* (const Float &b) {
			Su2	tmp;

			tmp.a[0] = a[0]*b;
			tmp.a[1] = a[1]*b;
			tmp.a[2] = a[2]*b;
			tmp.a[3] = a[3]*b;

			return	tmp;
		}

		Su2	operator*=(const Float &b) {
			(*this) = (*this)*b;
			return	(*this);
		}

		Float	Norm() {
			return	a[0]*a[0] + a[1]*a[1] + a[2]*a[2] + a[3]*a[3];
		}

		Float	Trace() {
			return	((Float) 2.)*a[0];
		}

		Float	SuperTrace() {	// Only different for vectorized datatypes
			return	((Float) 2.)*a[0];
		}

		Su2	Project() {
			return	(*this)/sqrt(Norm());
		}


		Su2&	GenHeat(const Float e2) {
			Su2   tmp;
			Float r, s, z;
			const Float oDet = sqrt(Norm());
			const Float uDet = 1./oDet;
			const Float e    = e2*oDet;
			const Float b    = exp(((Float) -2.) * e);

			do {
				r = Su2Rand::genRand();
				s = (((Float) 1.) - r) * b + r;
				tmp.a[0] = ((Float) 1.) + log(s) / e;

				r = Su2Rand::genRand();
				s = r*r;
				z = ((Float) 1.) - tmp.a[0]*tmp.a[0];
			}	while (s > z);

			s = sqrt(z);
			r = Su2Rand::genRand();

			tmp.a[3] = (((Float) 2.)*r - ((Float) 1.))*s;

			s = sqrt(std::abs(z - tmp.a[3]*tmp.a[3]));
			r = ((Float) 2.)*((Float) M_PI)*Su2Rand::genRand();

			tmp.a[1] = s*cos(r);
			tmp.a[2] = s*sin(r);

			a[0] *=  uDet;
			a[1] *= -uDet;
			a[2] *= -uDet;
			a[3] *= -uDet;

			(*this) = tmp*(*this);

			return	(*this);
		}

		Su2&	SetRandom(const Float eps = 1.) {
			#pragma unroll
			for (int i=0; i<4; i++) {
				Float r = Su2Rand::genRand();
				a[i] = ((Float) 2.*eps)*r - eps;
			}

			auto mod = sqrt(1./(a[0]*a[0] + a[1]*a[1] + a[2]*a[2] + a[3]*a[3]));

			a[0] *= mod;
			a[1] *= mod;
			a[2] *= mod;
			a[3] *= mod;

			return	(*this);
		}

		Su2&	Pert	(const Float eps) {
			Su2<Float> tmp;

			tmp.SetRandom(eps);

			(*this) *= tmp;
		}

		Su2	xPermute() {
			return	(*this);
		}

		Su2	yPermute() {
			return	(*this);
		}

		Su2	zPermute() {
			return	(*this);
		}

		Su2	tPermute() {
			return	(*this);
		}

		void	Print() {
			printf ("%+.3f %+.3f %+.3f %+.3f\n", a[0], a[1], a[2], a[3]);
		}

	};
#endif
