#ifndef	__AVX512F__
	#error("Building settings won't allow compilation for Avx-512")
#endif

#ifndef __SIMD
        #define __SIMD

	#define opCode_P(x,y,...) x ## _ ## y (__VA_ARGS__)
	#define opCode_N(x,y,...) opCode_P(x, y, __VA_ARGS__)
	#define opCode(x,...) opCode_N(_PREFIX_, x, __VA_ARGS__)

	#include <immintrin.h>
	#include "random/random.h"

	#define _MData_ __m512
	#define _MInt_  __m512i
	#define _MHnt_  __m256i
	#define _PREFIX_ _mm512
	#define _PREFXL_ _mm256
	#define opCodl(x,...) opCode_N(_PREFXL_, x, __VA_ARGS__)

        namespace Simd {

                constexpr size_t sAlign = 64;

		class	Simd_f;

		class	Mask_f {
			static constexpr size_t nData = sAlign/sizeof(int);

			private:

			__mmask16 data;

			public:

				Mask_f()			: data(0) {}
				Mask_f(const __mmask16 &in)	: data(in){}


			inline	Mask_f   operator& (const Mask_f &b) {
				return	opCode(kand, this->data, b.data);
			}

			inline	Mask_f  &operator&=(const Mask_f &b) {
				(*this) = (*this)&b;
				return	(*this);
			}

			inline	Mask_f   operator| (const Mask_f &b) {
				return	opCode(kor,  this->data, b.data);
			}

			inline	Mask_f  &operator|=(const Mask_f &b) {
				(*this) = (*this)|b;
				return  (*this);
			}

			inline	Mask_f   operator! () {
				return	opCode(kxor, this->data, 0b1111111111111111);
			}

			inline	int     Count() {
				int count = 0;

				for (int i=0; i<16; i++)
					count += (data>>i)&1;

				return  count;
			}

			void	Print(const char *str)  {
				printf("%s %x", str, this->data);
			}

			friend  class   Simd_f;
		};


		class	Simd_f	{
			private:

			_MData_	data;

			public:

                	static constexpr size_t	nData = sAlign/sizeof(float);
			static constexpr size_t	xWide = 2;
			static constexpr size_t	yWide = 2;
			static constexpr size_t	zWide = 2;
			static constexpr size_t	tWide = 2;

			typedef float sData;
			typedef Mask_f Mask;

				Simd_f() {
				data = opCode(setzero_ps);
			}

				Simd_f(float x0, float x1, float x2,  float x3,  float x4,  float x5,  float x6,  float x7,
				       float x8, float x9, float x10, float x11, float x12, float x13, float x14, float x15) {
				data = opCode(set_ps, x15, x14, x13, x12, x11, x10, x9, x8, x7, x6, x5, x4, x3, x2, x1, x0);
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

			inline	Simd_f	operator!() {
			        return  opCode(add_ps, opCode(permute_ps, this->data, 0b10110001), this->data);
			}

			inline	Simd_f	operator-() {
				return	opCode(sub_ps, opCode(setzero_ps), this->data);
			}

			inline	Simd_f	operator~() {
				return	opCode(mul_ps, this->data, opCode(set_ps, -1., 0., -1., 0., -1., 0., -1., 0., -1., 0., -1., 0., -1., 0., -1., 0.));
			}

			inline	Simd_f	operator^(const Mask_f &msk) {
				return	opCode(maskz_mov_ps, msk.data, this->data);
			}

			inline	Mask_f	operator>(const Simd_f &b) {
				return	opCode(cmp_ps_mask, this->data, b.data, _CMP_GT_OQ);
			}

			inline	Mask_f	operator>=(const Simd_f &b) {
				return	opCode(cmp_ps_mask, this->data, b.data, _CMP_GE_OQ);
			}

			inline	Mask_f	operator<(const Simd_f &b) {
				return	opCode(cmp_ps_mask, this->data, b.data, _CMP_LT_OQ);
			}

			inline	Mask_f	operator<=(const Simd_f &b) {
				return	opCode(cmp_ps_mask, this->data, b.data, _CMP_LE_OQ);
			}

			inline	Mask_f	operator==(const Simd_f &b) {
				return	opCode(cmp_ps_mask, this->data, b.data, _CMP_EQ_UQ);
			}

			inline	Simd_f	fma(const Simd_f &a, const Simd_f &b) {
				return	opCode(fmadd_ps, this->data, a.data, b.data);
			}

			inline	Simd_f	fms(const Simd_f &a, const Simd_f &b) {
				return	opCode(fmsub_ps, this->data, a.data, b.data);
			}

			inline	Simd_f	xPermute () {
				return	opCode(shuffle_f32x4, this->data, this->data, 0b01001110);
			}

			inline	Simd_f	yPermute () {
				return	opCode(shuffle_f32x4, this->data, this->data, 0b10110001);
			}

			inline	Simd_f	zPermute () {
				return	opCode(permute_ps, this->data, 0b01001110);
			}

			inline	Simd_f	tPermute () {
				return	opCode(permute_ps, this->data, 0b10110001);
			}

			inline	void	SetRandom () {
				(*this) = Simd_f(Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand(),
						 Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand(),
						 Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand(),
						 Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand(), Su2Rand::genRand());
			}

			inline	float	Sum () {
				return  opCode(reduce_add_ps, (*this).data);
			}

			inline	float&	operator[](int lane) {
				return	data[lane];
			}

			inline	_MData_&	raw() {
				return	data;
			}

			void	Print(const char *str) {
				printsVar(this->data, str);
			}

			friend	class	Mask_f;

			friend  Simd_f  sqrt    (const Simd_f&);
			friend  Simd_f  cos     (const Simd_f&);
			friend  Simd_f  sin     (const Simd_f&);
			friend  Simd_f  log     (const Simd_f&);
			friend  Simd_f  exp     (const Simd_f&);
		};

		Simd_f  sqrt    (const Simd_f&);
                Simd_f  cos     (const Simd_f&);
                Simd_f  sin     (const Simd_f&);
                Simd_f  log     (const Simd_f&);
                Simd_f  exp     (const Simd_f&);
	}
#endif
