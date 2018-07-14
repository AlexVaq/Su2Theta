#ifndef	__AVX512F__
	#error("Building settings won't allow compilation for Avx-512")
#endif

#ifndef __SIMD
        #define __SIMD

	#define opCode_P(x,y,...) x ## _ ## y (__VA_ARGS__)
	#define opCode_N(x,y,...) opCode_P(x, y, __VA_ARGS__)
	#define opCode(x,...) opCode_N(_PREFIX_, x, __VA_ARGS__)

	#include <immintrin.h>

	#define _MData_ __m512
	#define _MInt_  __m512i
	#define _MHnt_  __m256i
	#define _PREFIX_ _mm512
	#define _PREFXL_ _mm256
	#define opCodl(x,...) opCode_N(_PREFXL_, x, __VA_ARGS__)

        namespace Simd {

                constexpr size_t sAlign = 64;

		class	Simd_f	{
			private:

			_MData_	data;

			public:

				Simd_f(float x0, float x1, float x2, float x3, float x4, float x5, float x6, float x7, float x8, float x9, float x10, float x11, float x12, float x13, float x14, float x15) {
				data = opCode(set_ps, x15, x14, x13, x12, x11, x10, x9, x8, x7, x6, x5, x4, x3, x2, x1, x0);
			}

				Simd_f(float x) {
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
				opCode(store_ps, static_cast<void*>(memAddress));
			}

			void	stream (float *memAddress) {
				opCode(stream_ps, static_cast<void*>(memAddress), this->data);
			}

			inline	Simd_f&	operator+(Simd_f &x) {
				return	opCode(add_ps, this->data, x->data);
			}

			inline	Simd_f&	operator-(Simd_f &x) {
				return	opCode(sub_ps, this->data, x->data);
			}

			inline	Simd_f&	operator*(Simd_f &x) {
				return	opCode(mul_ps, this->data, x->data);
			}

			inline	Simd_f&	operator/(Simd_f &x) {
				return	opCode(div_ps, this->data, x->data);
			}

			inline	Simd_f&	operator+=(Simd_f &x) {
				return	(*this)+x;
			}

			inline	Simd_f&	operator-(Simd_f &x) {
				return	(*this)-x;
			}

			inline	Simd_f&	operator*(Simd_f &x) {
				return	(*this)*x;
			}

			inline	Simd_f&	operator/(Simd_f &x) {
				return	(*this)/x;
			}

			inline	Simd_f&	operator!() {
			        return  opCode(add_ps, opCode(permute_ps, this->data, 0b10110001), this->data);
			}

			inline	Simd_f&	operator~() {
				return	opCode(mul_ps, this->data, opCode(set_ps, -1., 0., -1., 0., -1., 0., -1., 0., -1., 0., -1., 0., -1., 0., -1., 0.));
			}

			inline	Simd_f&	fma(Simd_f &a, Simd_f &b) {
				return	opCode(fmadd_ps, this->data, a, b);
			}

			inline	Simd_f&	fms(Simd_f &a, Simd_f &b) {
				return	opCode(fmsub_ps, this->data, a, b);
			}

			inline	Simd_f&	xPermute () {
				return	opCode(shuffle_f32x4, this->data, this->data, 0b01001110);
			}

			inline	Simd_f&	yPermute () {
				return	opCode(shuffle_f32x4, this->data, this->data, 0b10110001);
			}

			inline	Simd_f&	zPermute () {
				return	opCode(permute_ps, this->data, 0b01001110);
			}

			inline	Simd_f&	tPermute () {
				return	opCode(permute_ps, this->data, 0b10110001);
			}

			inline	float	operator[](int lane) {
				return	data[lane];
			}
		}
	}
#endif
