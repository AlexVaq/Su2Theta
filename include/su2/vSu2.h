#ifndef	__VSU2_CLASS
	#define	__VSU2_CLASS

	#include "su2/su2.h"
	#include "simd/simd.h"
	#include "random/random.h"
	#include <algorithm>

	using namespace Simd;

	template<typename vFloat>
	class	vSu2	{
		private:

		vFloat	a[4];

		public:

		typedef		 vFloat        data;
		typedef	typename vFloat::sData sData;
		typedef	typename vFloat::Mask  Mask;
		typedef		 Su2<sData>    sClass;

		static constexpr size_t sWide = sizeof(vFloat)/sizeof(sData);
		static constexpr size_t xWide = vFloat::xWide;
		static constexpr size_t yWide = vFloat::yWide;
		static constexpr size_t zWide = vFloat::zWide;
		static constexpr size_t tWide = vFloat::tWide;

		vSu2	() 					     { }
		vSu2	(vFloat a0, vFloat a1, vFloat a2, vFloat a3) { a[0] = a0;   a[1] = a1;   a[2] = a2;   a[3] = a3; }
		vSu2	(sData  c)				     { a[0] = a[1] = a[2] = a[3] = vFloat(c); }
		vSu2	(sData * __restrict__ b)		     { a[0].Load(b); a[1].Load(b+sWide); a[2].Load(b+2*sWide); a[3].Load(b+3*sWide); }


		/*	Vectorized datatypes(CPU)	*/

		void	Save	(sData * __restrict__ out) {
			a[0].Save(out);
			a[1].Save(out+1*sWide);
			a[2].Save(out+2*sWide);
			a[3].Save(out+3*sWide);
		}

		void	SaveMask(Mask msk, sData * __restrict__ out) {
			a[0].SaveMask(msk, out);
			a[1].SaveMask(msk, out+1*sWide);
			a[2].SaveMask(msk, out+2*sWide);
			a[3].SaveMask(msk, out+3*sWide);
		}

		void	Stream	(sData * __restrict__ out) {
			a[0].Stream(out);
			a[1].Stream(out+1*sWide);
			a[2].Stream(out+2*sWide);
			a[3].Stream(out+3*sWide);
		}

		void	StreamMask(Mask msk, const vSu2 &in, sData * __restrict__ out) {
			a[0].StreamMask(msk, in.a[0], out);
			a[1].StreamMask(msk, in.a[1], out+1*sWide);
			a[2].StreamMask(msk, in.a[2], out+2*sWide);
			a[3].StreamMask(msk, in.a[3], out+3*sWide);
		}

		vFloat	operator% (const vSu2 &b) {
			return	a[0]*b.a[0] - a[1]*b.a[1] - a[2]*b.a[2] - a[3]*b.a[3];
		}

		vSu2	operator* (const vSu2 &b) {
			vSu2	tmp;

			tmp.a[0] = a[0]*b.a[0] - a[1]*b.a[1] - a[2]*b.a[2] - a[3]*b.a[3];
			tmp.a[1] = a[0]*b.a[1] + a[1]*b.a[0] - a[2]*b.a[3] + a[3]*b.a[2];
			tmp.a[2] = a[0]*b.a[2] + a[1]*b.a[3] + a[2]*b.a[0] - a[3]*b.a[1];
			tmp.a[3] = a[0]*b.a[3] - a[1]*b.a[2] + a[2]*b.a[1] + a[3]*b.a[0];

			return	tmp;
		}

		vSu2	operator*=(const vSu2 &b) {
			(*this) = (*this)*b;
			return	(*this);
		}

		vSu2	operator+ (const vSu2 &b) {
			vSu2	tmp;

			tmp.a[0] = a[0] + b.a[0];
			tmp.a[1] = a[1] + b.a[1];
			tmp.a[2] = a[2] + b.a[2];
			tmp.a[3] = a[3] + b.a[3];

			return	tmp;
		}

		vSu2	operator+=(const vSu2 &b) {
			(*this) = (*this)+b;
			return	(*this);
		}

		vSu2	operator- (const vSu2 &b) {
			vSu2	tmp;

			tmp.a[0] = a[0] - b.a[0];
			tmp.a[1] = a[1] - b.a[1];
			tmp.a[2] = a[2] - b.a[2];
			tmp.a[3] = a[3] - b.a[3];

			return	tmp;
		}

		vSu2	operator-=(const vSu2 &b) {
			(*this) = (*this)-b;
			return	(*this);
		}

		vSu2	operator! () {
			vSu2	tmp;

			tmp.a[0] =  a[0];
			tmp.a[1] = -a[1];
			tmp.a[2] = -a[2];
			tmp.a[3] = -a[3];

			return	tmp;
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
			(*this) = (*this)*b;
			return	(*this);
		}

		vSu2	operator* (const sData &b) {
			return	(*this)*vFloat(b);
		}

		vSu2	operator*=(const sData &b) {
			(*this) = (*this)*b;
			return	(*this);
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
			vSu2   tmp;
			vFloat r, s, z;
			vFloat oDet = sqrt(Norm());
			vFloat uDet = vFloat(1.)/oDet;
			vFloat e    = oDet*vFloat(e2);
			vFloat b    = exp(e*vFloat(-2.));
			vFloat aN;
			Mask   msk, tMsk;
			bool   notReady;

			do {
				r  = Su2Rand::vGenRand<vFloat>();//genVRand<vFloat>();//vFloat(0.5);
				s  = (vFloat(1.) - r) * b + r;
				aN = vFloat(1.) + log(s) / e;

				r  = Su2Rand::vGenRand<vFloat>();//genVRand<vFloat>();//vFloat(0.0);
				s  = r*r;
				z  = vFloat(1.) - aN*aN;

				msk = (s<=z);

				tmp.a[0] = (aN^msk) + (tmp.a[0]^(!msk));
				tMsk |= msk;
				notReady = (tMsk.Count() != vFloat::sWide) ? true : false;
			}	while (notReady);

			z    = vFloat(1.) - tmp.a[0]*tmp.a[0];
			s    = sqrt(z);
			r    = Su2Rand::vGenRand<vFloat>();//genVRand<vFloat>();//vFloat(0.5);

			tmp.a[3] = (vFloat(2.)*r - vFloat(1.))*s;

			s    = sqrt(abs(z - tmp.a[3]*tmp.a[3]));
			r    = vFloat(2.)*vFloat(M_PI)*Su2Rand::vGenRand<vFloat>();//genVRand<vFloat>();//(vFloat(0.5));

			tmp.a[1] = s*cos(r);
			tmp.a[2] = s*sin(r);

			a[0] *=  uDet;
			a[1] *= -uDet;
			a[2] *= -uDet;
			a[3] *= -uDet;

			(*this) = tmp*(*this);

			return	(*this);
		}

		vSu2	Reflect	(vSu2 &&staple) {
			auto sNrm = staple.Norm()*0.5;
			auto lmba = ((*this)%staple)/sNrm;
			return	((!staple)*lmba) - (*this);
		}

		vSu2&	SetRandom(const sData eps = 1.) {
			#pragma unroll
			for (int i=0; i<4; i++) {
				vFloat r = Su2Rand::vGenRand<vFloat>();//genVRand<vFloat>();
				a[i] = vFloat(2.*eps)*r - vFloat(eps);
			}

			auto mod = sqrt(vFloat(1.)/(a[0]*a[0] + a[1]*a[1] + a[2]*a[2] + a[3]*a[3]));

			a[0] *= mod;
			a[1] *= mod;
			a[2] *= mod;
			a[3] *= mod;

			return	(*this);
		}

		vSu2	Pert	(const typename vFloat::sData eps) {
			vSu2<vFloat> tmp;

			tmp.SetRandom(eps);

			return	tmp*(*this);
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

		Su2<sData>	Extract(const int lane) {
			auto a0 = a[0][lane];
			auto a1 = a[1][lane];
			auto a2 = a[2][lane];
			auto a3 = a[3][lane];

			return	Su2<sData>(a0, a1, a2, a3);
		}

		void	Print() {
			a[0].Print("0: ");
			a[1].Print("1: ");
			a[2].Print("2: ");
			a[3].Print("3: ");
		}
	};
#endif
