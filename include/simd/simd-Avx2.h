#if	!defined(__AVX2__) && !defined(__AVX__)
	#error("Building settings won't allow compilation for Avx/Avx-2")
#endif

#ifndef __SIMD
        #define __SIMD

	#define opCode_P(x,y,...) x ## _ ## y (__VA_ARGS__)
	#define opCode_N(x,y,...) opCode_P(x, y, __VA_ARGS__)
	#define opCode(x,...) opCode_N(_PREFIX_, x, __VA_ARGS__)

	#include <immintrin.h>
	#include "random/random.h"

	#define _MData_ __m256
	#define _MInt_  __m256i
	#define _MHnt_  __m128i
	#define _PREFIX_ _mm256
	#define _PREFXL_ _mm128
	#define opCodl(x,...) opCode_N(_PREFXL_, x, __VA_ARGS__)

        namespace Simd {

                constexpr size_t sAlign = 32;

		class	Simd_f	{
			private:

			_MData_	data;

			public:

                	static constexpr size_t nData = sAlign/sizeof(float);
			typedef float sData;

				Simd_f() {
				data = opCode(setzero_ps);
			}

				Simd_f(float x0, float x1, float x2, float x3, float x4, float x5, float x6, float x7) {
				data = opCode(set_ps, x7, x6, x5, x4, x3, x2, x1, x0);
			}

				Simd_f(float x0) {
				data = opCode(set1_ps, x0);
			}

				Simd_f(float *memAddress) {
				data = opCode(load_ps, memAddress);
			}

				Simd_f(const _MData_ &in) : data(in) {};
				Simd_f(_MData_ &&in) : data(std::move(in)) {};

			inline	Simd_f&	operator=(const _MData_ &in) {
				data = in;
			}

			inline	Simd_f&	operator=(_MData_ &&in) {
				data = std::move(in);
			}

			void	save  (float *memAddress) {
				opCode(store_ps,  memAddress, this->data);
			}

			void	stream (float *memAddress) {
				opCode(stream_ps, memAddress, this->data);
			}

			inline	Simd_f	operator+(const Simd_f &x) {
				return	opCode(add_ps, this->data, x.data);
			}

			inline	Simd_f	operator-(const Simd_f &x) {
				return	opCode(sub_ps, this->data, x.data);
			}

			inline	Simd_f	operator*(const Simd_f &x) {
				return	opCode(mul_ps, this->data, x.data);
			}

			inline	Simd_f	operator/(const Simd_f &x) {
				return	opCode(div_ps, this->data, x.data);
			}

			inline	Simd_f &operator+=(const Simd_f &x) {
				(*this) = (*this)+x;
				return	(*this);
			}

			inline	Simd_f &operator-=(const Simd_f &x) {
				(*this) = (*this)-x;
				return	(*this);
			}

			inline	Simd_f &operator*=(const Simd_f &x) {
				(*this) = (*this)*x;
				return	(*this);
			}

			inline	Simd_f &operator/=(const Simd_f &x) {
				(*this) = (*this)/x;
				return	(*this);
			}

			inline	Simd_f	operator-() {
				return	opCode(sub_ps, opCode(setzero_ps), this->data);
			}

			inline	Simd_f	operator!() {
				return	opCode(add_ps, opCode(permute_ps, this->data, 0b10110001), this->data);
			}

			inline	Simd_f	operator~() {
				return	opCode(mul_ps, this->data, opCode(set_ps, -1., 0., -1., 0., -1., 0., -1., 0.));
			}

			inline	Simd_f	fma(const Simd_f &a, const Simd_f &b) {
				return	opCode(fmadd_ps, this->data, a.data, b.data);
			}

			inline	Simd_f	fms(const Simd_f &a, const Simd_f &b) {
				return	opCode(fmsub_ps, this->data, a.data, b.data);
			}

			inline	Simd_f	xPermute () {
				return	(*this);
			}

			inline	Simd_f	yPermute () {
				return	opCode(permute2f128_ps, this->data, this->data, 0b00000001);
			}

			inline	Simd_f	zPermute () {
				return	opCode(permute_ps, this->data, 0b01001110);
			}

			inline	Simd_f	tPermute () {
				return	opCode(permute_ps, this->data, 0b10110001);
			}

			inline	void	SetRandom () {
				(*this) = Simd_f(Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand(),
						 Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand());
			}

			inline	float	Sum () {
				return	opCode(hadd_ps, opCode(hadd_ps,
						opCode(hadd_ps, (*this).data, opCode(permute2f128_ps, (*this).data, (*this).data, 0b00000001)), (*this).data), (*this).data)[0];
			}

			inline	float&	operator[](int lane) {
				return	data[lane];
			}

			inline	_MData_&	raw() {
				return	data;
			}

			friend  Simd_f  sqrt    (const Simd_f&);
			friend  Simd_f  cos     (const Simd_f&);
			friend  Simd_f  sin     (const Simd_f&);
		};
	}
#endif
