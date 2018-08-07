#ifndef	__TABLEINTRINSICS
#define	__TABLEINTRINSICS

#include<cmath>
#include<array>

#define opCode_P(x,y,...) x ## _ ## y (__VA_ARGS__)
#define opCode_N(x,y,...) opCode_P(x, y, __VA_ARGS__)
#define opCode(x,...) opCode_N(_PREFIX_, x, __VA_ARGS__)

#include <immintrin.h>

#ifdef	__AVX512F__
	#define _MData_ __m512d
	#define	_MInt_  __m512i
	#define	_MHnt_  __m256i
#elif   defined(__AVX__)
	#define _MData_ __m256d
	#define	_MInt_  __m256i
	#define	_MHnt_  __m128i
#else
	#define _MData_ __m128d
	#define	_MInt_  __m128i
#endif

#if	defined(__AVX512F__)
	#define	_PREFIX_ _mm512
	#define	_PREFXL_ _mm256
	#define opCodl(x,...) opCode_N(_PREFXL_, x, __VA_ARGS__)
#else
	#if not defined(__AVX__) and not defined(__AVX2__)
		#define	_PREFIX_ _mm
	#else
		#define	_PREFIX_ _mm256
		#define	_PREFXL_ _mm
		#define opCodl(x,...) opCode_N(_PREFXL_, x, __VA_ARGS__)
	#endif
#endif

#define	M_PI2	(M_PI *M_PI)
#define	M_PI4	(M_PI2*M_PI2)
#define	M_PI6	(M_PI4*M_PI2)
/* Only single precision for the moment */

#ifndef	__INTEL_COMPILER

constexpr double Inf_d = __builtin_inf();
constexpr double Nan_d = __builtin_nan("");//0xFFFFF");

constexpr size_t m32Mk = 0b0000000000000000000000000001111100000000000000000000000000011111;

constexpr double PiA_d = -3.1415926218032836914;
constexpr double PiB_d = -3.1786509424591713469e-08;
constexpr double PiC_d = -1.2246467864107188502e-16;
constexpr double PiD_d = -1.2736634327021899816e-24;

constexpr double s0_d  = -7.97255955009037868891952e-18;
constexpr double s1_d  =  2.81009972710863200091251e-15;
constexpr double s2_d  = -7.64712219118158833288484e-13;
constexpr double s3_d  =  1.60590430605664501629054e-10;
constexpr double s4_d  = -2.50521083763502045810755e-08;
constexpr double s5_d  =  2.75573192239198747630416e-06;
constexpr double s6_d  = -0.000198412698412696162806809;
constexpr double s7_d  =  0.00833333333333332974823815;
constexpr double s8_d  = -0.166666666666666657414808;
#ifdef	__AVX512F__
constexpr _MData_ rPid	    = { M_1_PI/(1<<24), M_1_PI/(1<<24), M_1_PI/(1<<24), M_1_PI/(1<<24), M_1_PI/(1<<24), M_1_PI/(1<<24), M_1_PI/(1<<24), M_1_PI/(1<<24) };
constexpr _MData_ dPid	    = { M_1_PI/(1<<23), M_1_PI/(1<<23), M_1_PI/(1<<23), M_1_PI/(1<<23), M_1_PI/(1<<23), M_1_PI/(1<<23), M_1_PI/(1<<23), M_1_PI/(1<<23) };
constexpr _MData_ rCte      = {   16777216.,   16777216.,   16777216.,   16777216.,   16777216.,   16777216.,   16777216.,   16777216. };
constexpr _MData_ oPid      = {     1./M_PI,     1./M_PI,     1./M_PI,     1./M_PI,     1./M_PI,     1./M_PI,     1./M_PI,     1./M_PI };
constexpr _MData_ zeroNegd  = {        -0.0,        -0.0,        -0.0,        -0.0,        -0.0,        -0.0,        -0.0,        -0.0 };
constexpr _MData_ dHlf      = {         0.5,         0.5,         0.5,         0.5,         0.5,         0.5,         0.5,         0.5 };
constexpr _MData_ dOne      = {         1.0,         1.0,         1.0,         1.0,         1.0,         1.0,         1.0,         1.0 };
constexpr _MData_ TriMaxd   = {        1e15,        1e15,        1e15,        1e15,        1e15,        1e15,        1e15,        1e15 };
constexpr _MData_ dInf      = {       Inf_d,       Inf_d,       Inf_d,       Inf_d,       Inf_d,       Inf_d,       Inf_d,       Inf_d };
constexpr _MData_ dNan      = {       Nan_d,       Nan_d,       Nan_d,       Nan_d,       Nan_d,       Nan_d,       Nan_d,       Nan_d };
constexpr _MData_ PiAd      = {       PiA_d,       PiA_d,       PiA_d,       PiA_d,       PiA_d,       PiA_d,       PiA_d,       PiA_d };
constexpr _MData_ PiBd      = {       PiB_d,       PiB_d,       PiB_d,       PiB_d,       PiB_d,       PiB_d,       PiB_d,       PiB_d };
constexpr _MData_ PiCd      = {       PiC_d,       PiC_d,       PiC_d,       PiC_d,       PiC_d,       PiC_d,       PiC_d,       PiC_d };
constexpr _MData_ PiDd      = {       PiD_d,       PiD_d,       PiD_d,       PiD_d,       PiD_d,       PiD_d,       PiD_d,       PiD_d };
constexpr _MData_ hPiAd     = {   0.5*PiA_d,   0.5*PiA_d,   0.5*PiA_d,   0.5*PiA_d,   0.5*PiA_d,   0.5*PiA_d,   0.5*PiA_d,   0.5*PiA_d };
constexpr _MData_ hPiBd     = {   0.5*PiB_d,   0.5*PiB_d,   0.5*PiB_d,   0.5*PiB_d,   0.5*PiB_d,   0.5*PiB_d,   0.5*PiB_d,   0.5*PiB_d };
constexpr _MData_ hPiCd     = {   0.5*PiC_d,   0.5*PiC_d,   0.5*PiC_d,   0.5*PiC_d,   0.5*PiC_d,   0.5*PiC_d,   0.5*PiC_d,   0.5*PiC_d };
constexpr _MData_ hPiDd     = {   0.5*PiD_d,   0.5*PiD_d,   0.5*PiD_d,   0.5*PiD_d,   0.5*PiD_d,   0.5*PiD_d,   0.5*PiD_d,   0.5*PiD_d };
constexpr _MData_ s0d       = {        s0_d,        s0_d,        s0_d,        s0_d,        s0_d,        s0_d,        s0_d,        s0_d };
constexpr _MData_ s1d       = {        s1_d,        s1_d,        s1_d,        s1_d,        s1_d,        s1_d,        s1_d,        s1_d };
constexpr _MData_ s2d       = {        s2_d,        s2_d,        s2_d,        s2_d,        s2_d,        s2_d,        s2_d,        s2_d };
constexpr _MData_ s3d       = {        s3_d,        s3_d,        s3_d,        s3_d,        s3_d,        s3_d,        s3_d,        s3_d };
constexpr _MData_ s4d       = {        s4_d,        s4_d,        s4_d,        s4_d,        s4_d,        s4_d,        s4_d,        s4_d };
constexpr _MData_ s5d       = {        s5_d,        s5_d,        s5_d,        s5_d,        s5_d,        s5_d,        s5_d,        s5_d };
constexpr _MData_ s6d       = {        s6_d,        s6_d,        s6_d,        s6_d,        s6_d,        s6_d,        s6_d,        s6_d };
constexpr _MData_ s7d       = {        s7_d,        s7_d,        s7_d,        s7_d,        s7_d,        s7_d,        s7_d,        s7_d };
constexpr _MData_ s8d       = {        s8_d,        s8_d,        s8_d,        s8_d,        s8_d,        s8_d,        s8_d,        s8_d };
constexpr _MInt_  iZero     = {           0,           0,           0,           0,           0,           0,           0,           0 };
constexpr _MInt_  one       = {  4294967297,  4294967297,  4294967297,  4294967297,  4294967297,  4294967297,  4294967297,  4294967297 };
constexpr _MInt_  two       = {  8589934594,  8589934594,  8589934594,  8589934594,  8589934594,  8589934594,  8589934594,  8589934594 };
constexpr _MHnt_  hOne      = {  4294967297,  4294967297,  4294967297,  4294967297 };
constexpr _MHnt_  hTwo      = {  8589934594,  8589934594,  8589934594,  8589934594 };
constexpr _MHnt_  iZerh     = {           0,           0,           0,           0 };
constexpr _MInt_  m32Mask   = {       m32Mk,       m32Mk,       m32Mk,       m32Mk,       m32Mk,       m32Mk,       m32Mk,       m32Mk };
#elif   defined(__AVX__)
constexpr _MData_ rPid	    = { M_1_PI/(1<<24), M_1_PI/(1<<24), M_1_PI/(1<<24), M_1_PI/(1<<24) };
constexpr _MData_ dPid	    = { M_1_PI/(1<<23), M_1_PI/(1<<23), M_1_PI/(1<<23), M_1_PI/(1<<23) };
constexpr _MData_ rCte      = {     (1<<24),     (1<<24),     (1<<24),     (1<<24) };
constexpr _MData_ dCte      = {     (1<<23),     (1<<23),     (1<<23),     (1<<23) };
constexpr _MData_ oPid      = {      M_1_PI,      M_1_PI,      M_1_PI,      M_1_PI };
constexpr _MData_ zeroNegd  = {        -0.0,        -0.0,        -0.0,        -0.0 };
constexpr _MData_ dHlf      = {         0.5,         0.5,         0.5,         0.5 };
constexpr _MData_ dOne      = {         1.0,         1.0,         1.0,         1.0 };
#ifdef	__FMA__
constexpr _MData_ TriMaxd   = {        1e15,        1e15,        1e15,        1e15 };
#else
constexpr _MData_ TriMaxd   = {        1e12,        1e12,        1e12,        1e12 };
#endif
constexpr _MData_ dInf      = {       Inf_d,       Inf_d,       Inf_d,       Inf_d };
constexpr _MData_ dNan      = {       Nan_d,       Nan_d,       Nan_d,       Nan_d };
constexpr _MData_ PiAd      = {       PiA_d,       PiA_d,       PiA_d,       PiA_d };
constexpr _MData_ PiBd      = {       PiB_d,       PiB_d,       PiB_d,       PiB_d };
constexpr _MData_ PiCd      = {       PiC_d,       PiC_d,       PiC_d,       PiC_d };
constexpr _MData_ PiDd      = {       PiD_d,       PiD_d,       PiD_d,       PiD_d };
constexpr _MData_ hPiAd     = {   0.5*PiA_d,   0.5*PiA_d,   0.5*PiA_d,   0.5*PiA_d };
constexpr _MData_ hPiBd     = {   0.5*PiB_d,   0.5*PiB_d,   0.5*PiB_d,   0.5*PiB_d };
constexpr _MData_ hPiCd     = {   0.5*PiC_d,   0.5*PiC_d,   0.5*PiC_d,   0.5*PiC_d };
constexpr _MData_ hPiDd     = {   0.5*PiD_d,   0.5*PiD_d,   0.5*PiD_d,   0.5*PiD_d };
constexpr _MData_ s0d       = {        s0_d,        s0_d,        s0_d,        s0_d };
constexpr _MData_ s1d       = {        s1_d,        s1_d,        s1_d,        s1_d };
constexpr _MData_ s2d       = {        s2_d,        s2_d,        s2_d,        s2_d };
constexpr _MData_ s3d       = {        s3_d,        s3_d,        s3_d,        s3_d };
constexpr _MData_ s4d       = {        s4_d,        s4_d,        s4_d,        s4_d };
constexpr _MData_ s5d       = {        s5_d,        s5_d,        s5_d,        s5_d };
constexpr _MData_ s6d       = {        s6_d,        s6_d,        s6_d,        s6_d };
constexpr _MData_ s7d       = {        s7_d,        s7_d,        s7_d,        s7_d };
constexpr _MData_ s8d       = {        s8_d,        s8_d,        s8_d,        s8_d };
constexpr _MInt_  iZero     = {           0,           0,           0,           0 };
constexpr _MInt_  one       = {  4294967297,  4294967297,  4294967297,  4294967297 };
constexpr _MInt_  two       = {  8589934594,  8589934594,  8589934594,  8589934594 };
constexpr _MHnt_  iZerh     = {           0,           0 };
constexpr _MHnt_  hOne      = {  4294967297,  4294967297 };
constexpr _MHnt_  hTwo      = {  8589934594,  8589934594 };
constexpr _MInt_  m32Mask   = {       m32Mk,       m32Mk,       m32Mk,       m32Mk };
#else
constexpr _MData_ rPid	    = { M_1_PI/(1<<24), M_1_PI/(1<<24) };
constexpr _MData_ dPid	    = { M_1_PI/(1<<23), M_1_PI/(1<<23) };
constexpr _MData_ rCte      = {   (1 << 24),   (1 << 24) };
constexpr _MData_ dCte      = {   (1 << 23),   (1 << 23) };
constexpr _MData_ oPid      = {      M_1_PI,      M_1_PI };
constexpr _MData_ zeroNegd  = {        -0.0,        -0.0 };
constexpr _MData_ dHlf      = {         0.5,         0.5 };
constexpr _MData_ dOne      = {         1.0,         1.0 };
constexpr _MData_ TriMaxd   = {        1e15,        1e15 };
constexpr _MData_ dInf      = {       Inf_d,       Inf_d };
constexpr _MData_ dNan      = {       Nan_d,       Nan_d };
constexpr _MData_ PiAd      = {       PiA_d,       PiA_d };
constexpr _MData_ PiBd      = {       PiB_d,       PiB_d };
constexpr _MData_ PiCd      = {       PiC_d,       PiC_d };
constexpr _MData_ PiDd      = {       PiD_d,       PiD_d };
constexpr _MData_ hPiAd     = {   0.5*PiA_d,   0.5*PiA_d };
constexpr _MData_ hPiBd     = {   0.5*PiB_d,   0.5*PiB_d };
constexpr _MData_ hPiCd     = {   0.5*PiC_d,   0.5*PiC_d };
constexpr _MData_ hPiDd     = {   0.5*PiD_d,   0.5*PiD_d };
constexpr _MData_ s0d       = {        s0_d,        s0_d };
constexpr _MData_ s1d       = {        s1_d,        s1_d };
constexpr _MData_ s2d       = {        s2_d,        s2_d };
constexpr _MData_ s3d       = {        s3_d,        s3_d };
constexpr _MData_ s4d       = {        s4_d,        s4_d };
constexpr _MData_ s5d       = {        s5_d,        s5_d };
constexpr _MData_ s6d       = {        s6_d,        s6_d };
constexpr _MData_ s7d       = {        s7_d,        s7_d };
constexpr _MData_ s8d       = {        s8_d,        s8_d };
constexpr _MInt_  iZero     = {           0,           0 };
constexpr _MInt_  one       = {  4294967297,  4294967297 };
constexpr _MInt_  two       = {  8589934594,  8589934594 };
constexpr _MInt_  m32Mask   = {       m32Mk,       m32Mk };
#endif

#ifdef	__AVX__
inline void printhVar(_MHnt_ d, const char *name)
#else
inline void printhVar(_MInt_ d, const char *name)
#endif
{
	printf ("%s", name);
#if	defined(__AVX512F__)
	int r[8] __attribute((aligned(32)));
	opCodl(store_si256, ((_MHnt_ *)r), d);
	for (int i=0; i<8; i++)
#elif	defined(__AVX__)
	int r[4] __attribute((aligned(16)));
	opCodl(store_si128, ((_MHnt_ *)r), d);
	for (int i=0; i<4; i++)
#else
	int r[4] __attribute((aligned(16)));
	opCode(store_si128, ((_MInt_ *)r), d);
	for (int i=0; i<4; i++)
#endif
		printf(" %d", r[i]);
	printf("\n");
}

inline void printdVar(_MData_ d, const char *name) {

	printf ("%s", name);
#if	defined(__AVX512F__)
	for (int i=0; i<8; i++)
#elif	defined(__AVX__)
	for (int i=0; i<4; i++)
#else
	for (int i=0; i<2; i++)
#endif
		printf(" %+lf", d[i]);
	printf("\n");
}

#undef	_MData_

#if	defined(__AVX512F__)
	#define	_MData_ __m512
#elif	defined(__AVX__)
	#define	_MData_ __m256
#else
	#define	_MData_ __m128
#endif 
#endif
#ifndef	__INTEL_COMPILER

constexpr float Inf_f = __builtin_inff();
constexpr float Nan_f = __builtin_nanf("0x3FFFFF");

/*	For the exponential	*/
constexpr float th1_f = 220.4208f;
constexpr float th2_f = 2.9802322E-8f;
constexpr float ivL_f = 46.16624f;
constexpr float L1_f  = 0.021660805f;
constexpr float L2_f  = 4.464396e-8f;
constexpr float A1_f  = 0.50000405f;
constexpr float A2_f  = 0.16666764f;

constexpr float Sl00_f = 1.0000000f;
constexpr float Sl01_f = 1.0218964f;
constexpr float Sl02_f = 1.0442734f;
constexpr float Sl03_f = 1.0671387f;
constexpr float Sl04_f = 1.0905075f;
constexpr float Sl05_f = 1.1143799f;
constexpr float Sl06_f = 1.1387863f;
constexpr float Sl07_f = 1.1637192f;
constexpr float Sl08_f = 1.1892014f;
constexpr float Sl09_f = 1.2152405f;
constexpr float Sl10_f = 1.2418518f;
constexpr float Sl11_f = 1.2690506f;
constexpr float Sl12_f = 1.2968369f;
constexpr float Sl13_f = 1.3252335f;
constexpr float Sl14_f = 1.3542480f;
constexpr float Sl15_f = 1.3839035f;
constexpr float Sl16_f = 1.4142075f;
constexpr float Sl17_f = 1.4451752f;
constexpr float Sl18_f = 1.4768219f;
constexpr float Sl19_f = 1.5091629f;
constexpr float Sl20_f = 1.5422058f;
constexpr float Sl21_f = 1.5759735f;
constexpr float Sl22_f = 1.6104889f;
constexpr float Sl23_f = 1.6457520f;
constexpr float Sl24_f = 1.6817856f;
constexpr float Sl25_f = 1.7186127f;
constexpr float Sl26_f = 1.7562485f;
constexpr float Sl27_f = 1.7947083f;
constexpr float Sl28_f = 1.8340073f;
constexpr float Sl29_f = 1.8741608f;
constexpr float Sl30_f = 1.9151993f;
constexpr float Sl31_f = 1.9571381f;

constexpr float St00_f = 0.0000000f;
constexpr float St01_f = 7.8634940e-07f;
constexpr float St02_f = 4.0596257e-07f;
constexpr float St03_f = 1.7288019e-06f;
constexpr float St04_f = 2.2534104e-07f;
constexpr float St05_f = 6.8597833e-06f;
constexpr float St06_f = 2.3188388e-06f;
constexpr float St07_f = 5.6815315e-06f;
constexpr float St08_f = 5.7600223e-06f;
constexpr float St09_f = 6.8814647e-06f;
constexpr float St10_f = 6.0054331e-06f;
constexpr float St11_f = 3.5904719e-07f;
constexpr float St12_f = 2.7016238e-06f;
constexpr float St13_f = 3.1836871e-06f;
constexpr float St14_f = 7.5000621e-06f;
constexpr float St15_f = 6.3785460e-06f;
constexpr float St16_f = 6.1038768e-06f;
constexpr float St17_f = 5.6360786e-06f;
constexpr float St18_f = 4.2465254e-06f;
constexpr float St19_f = 1.5247614e-06f;
constexpr float St20_f = 5.0148610e-06f;
constexpr float St21_f = 7.3343658e-06f;
constexpr float St22_f = 1.4403477e-06f;
constexpr float St23_f = 3.5250289e-06f;
constexpr float St24_f = 7.2470111e-06f;
constexpr float St25_f = 6.6272241e-06f;
constexpr float St26_f = 3.6862523e-06f;
constexpr float St27_f = 8.2304996e-07f;
constexpr float St28_f = 8.2322578e-07f;
constexpr float St29_f = 6.8675085e-06f;
constexpr float St30_f = 7.2816119e-06f;
constexpr float St31_f = 6.0626521e-06f;

/*	For the trigonometric functions	*/

constexpr float PiA_f = -3.140625f;
constexpr float PiB_f = -0.0009670257568359375f;
constexpr float PiC_f = -6.2771141529083251953e-07f;
constexpr float PiD_f = -1.2154201256553420762e-10f;

constexpr float s0_f  =  2.6083159809786593541503e-06f;
constexpr float s1_f  = -0.0001981069071916863322258f;
constexpr float s2_f  =  0.00833307858556509017944336f;
constexpr float s3_f  = -0.166666597127914428710938f;

#ifdef	__AVX512F__
constexpr _MData_ oPif      = {     1./M_PI,     1./M_PI,     1./M_PI,     1./M_PI,     1./M_PI,     1./M_PI,     1./M_PI,     1./M_PI,
				    1./M_PI,     1./M_PI,     1./M_PI,     1./M_PI,     1./M_PI,     1./M_PI,     1./M_PI,     1./M_PI };
constexpr _MData_ zeroNegf  = {       -0.0f,       -0.0f,       -0.0f,       -0.0f,       -0.0f,       -0.0f,       -0.0f,       -0.0f,
				      -0.0f,       -0.0f,       -0.0f,       -0.0f,       -0.0f,       -0.0f,       -0.0f,       -0.0f };
constexpr _MData_ fHlf      = {        0.5f,        0.5f,        0.5f,        0.5f,        0.5f,        0.5f,        0.5f,        0.5f,
				       0.5f,        0.5f,        0.5f,        0.5f,        0.5f,        0.5f,        0.5f,        0.5f };
constexpr _MData_ TriMaxf   = {         1e7,         1e7,         1e7,         1e7,         1e7,         1e7,         1e7,         1e7,
					1e7,         1e7,         1e7,         1e7,         1e7,         1e7,         1e7,         1e7 };
constexpr _MData_ fInf      = {       Inf_f,       Inf_f,       Inf_f,       Inf_f,       Inf_f,       Inf_f,       Inf_f,       Inf_f,
				      Inf_f,       Inf_f,       Inf_f,       Inf_f,       Inf_f,       Inf_f,       Inf_f,       Inf_f };
constexpr _MData_ fNan      = {       Nan_f,       Nan_f,       Nan_f,       Nan_f,       Nan_f,       Nan_f,       Nan_f,       Nan_f,
				      Nan_f,       Nan_f,       Nan_f,       Nan_f,       Nan_f,       Nan_f,       Nan_f,       Nan_f };
constexpr _MData_ PiAf      = {       PiA_f,       PiA_f,       PiA_f,       PiA_f,       PiA_f,       PiA_f,       PiA_f,       PiA_f,
				      PiA_f,       PiA_f,       PiA_f,       PiA_f,       PiA_f,       PiA_f,       PiA_f,       PiA_f };
constexpr _MData_ PiBf      = {       PiB_f,       PiB_f,       PiB_f,       PiB_f,       PiB_f,       PiB_f,       PiB_f,       PiB_f,
				      PiB_f,       PiB_f,       PiB_f,       PiB_f,       PiB_f,       PiB_f,       PiB_f,       PiB_f };
constexpr _MData_ PiCf      = {       PiC_f,       PiC_f,       PiC_f,       PiC_f,       PiC_f,       PiC_f,       PiC_f,       PiC_f,
				      PiC_f,       PiC_f,       PiC_f,       PiC_f,       PiC_f,       PiC_f,       PiC_f,       PiC_f };
constexpr _MData_ PiDf      = {       PiD_f,       PiD_f,       PiD_f,       PiD_f,       PiD_f,       PiD_f,       PiD_f,       PiD_f,
				      PiD_f,       PiD_f,       PiD_f,       PiD_f,       PiD_f,       PiD_f,       PiD_f,       PiD_f };
constexpr _MData_ hPiAf     = {  0.5f*PiA_f,  0.5f*PiA_f,  0.5f*PiA_f,  0.5f*PiA_f,  0.5f*PiA_f,  0.5f*PiA_f,  0.5f*PiA_f,  0.5f*PiA_f,
				 0.5f*PiA_f,  0.5f*PiA_f,  0.5f*PiA_f,  0.5f*PiA_f,  0.5f*PiA_f,  0.5f*PiA_f,  0.5f*PiA_f,  0.5f*PiA_f };
constexpr _MData_ hPiBf     = {  0.5f*PiB_f,  0.5f*PiB_f,  0.5f*PiB_f,  0.5f*PiB_f,  0.5f*PiB_f,  0.5f*PiB_f,  0.5f*PiB_f,  0.5f*PiB_f,
				 0.5f*PiB_f,  0.5f*PiB_f,  0.5f*PiB_f,  0.5f*PiB_f,  0.5f*PiB_f,  0.5f*PiB_f,  0.5f*PiB_f,  0.5f*PiB_f };
constexpr _MData_ hPiCf     = {  0.5f*PiC_f,  0.5f*PiC_f,  0.5f*PiC_f,  0.5f*PiC_f,  0.5f*PiC_f,  0.5f*PiC_f,  0.5f*PiC_f,  0.5f*PiC_f,
				 0.5f*PiC_f,  0.5f*PiC_f,  0.5f*PiC_f,  0.5f*PiC_f,  0.5f*PiC_f,  0.5f*PiC_f,  0.5f*PiC_f,  0.5f*PiC_f };
constexpr _MData_ hPiDf     = {  0.5f*PiD_f,  0.5f*PiD_f,  0.5f*PiD_f,  0.5f*PiD_f,  0.5f*PiD_f,  0.5f*PiD_f,  0.5f*PiD_f,  0.5f*PiD_f,
				 0.5f*PiD_f,  0.5f*PiD_f,  0.5f*PiD_f,  0.5f*PiD_f,  0.5f*PiD_f,  0.5f*PiD_f,  0.5f*PiD_f,  0.5f*PiD_f };
constexpr _MData_ s0f       = {        s0_f,        s0_f,        s0_f,        s0_f,        s0_f,        s0_f,        s0_f,        s0_f,
				       s0_f,        s0_f,        s0_f,        s0_f,        s0_f,        s0_f,        s0_f,        s0_f };
constexpr _MData_ s1f       = {        s1_f,        s1_f,        s1_f,        s1_f,        s1_f,        s1_f,        s1_f,        s1_f,
				       s1_f,        s1_f,        s1_f,        s1_f,        s1_f,        s1_f,        s1_f,        s1_f };
constexpr _MData_ s2f       = {        s2_f,        s2_f,        s2_f,        s2_f,        s2_f,        s2_f,        s2_f,        s2_f,
				       s2_f,        s2_f,        s2_f,        s2_f,        s2_f,        s2_f,        s2_f,        s2_f };
constexpr _MData_ s3f       = {        s3_f,        s3_f,        s3_f,        s3_f,        s3_f,        s3_f,        s3_f,        s3_f,
				       s3_f,        s3_f,        s3_f,        s3_f,        s3_f,        s3_f,        s3_f,        s3_f };

constexpr _MData_ vTh1f     = {       th1_f,       th1_f,       th1_f,       th1_f,       th1_f,       th1_f,       th1_f,       th1_f,
				      th1_f,       th1_f,       th1_f,       th1_f,       th1_f,       th1_f,       th1_f,       th1_f };
constexpr _MData_ vTh2f     = {       th2_f,       th2_f,       th2_f,       th2_f,       th2_f,       th2_f,       th2_f,       th2_f,
				      th2_f,       th2_f,       th2_f,       th2_f,       th2_f,       th2_f,       th2_f,       th2_f };
constexpr _MData_ vL1f      = {        L1_f,        L1_f,        L1_f,        L1_f,        L1_f,        L1_f,        L1_f,        L1_f,
			 	       L1_f,        L1_f,        L1_f,        L1_f,        L1_f,        L1_f,        L1_f,        L1_f };
constexpr _MData_ vL2f      = {       -L2_f,       -L2_f,       -L2_f,       -L2_f,       -L2_f,       -L2_f,       -L2_f,       -L2_f,
				      -L2_f,       -L2_f,       -L2_f,       -L2_f,       -L2_f,       -L2_f,       -L2_f,       -L2_f };
constexpr _MData_ vA1f      = {        A1_f,        A1_f,        A1_f,        A1_f,        A1_f,        A1_f,        A1_f,        A1_f,
			 	       A1_f,        A1_f,        A1_f,        A1_f,        A1_f,        A1_f,        A1_f,        A1_f };
constexpr _MData_ vA2f      = {        A2_f,        A2_f,        A2_f,        A2_f,        A2_f,        A2_f,        A2_f,        A2_f,
				       A2_f,        A2_f,        A2_f,        A2_f,        A2_f,        A2_f,        A2_f,        A2_f };
constexpr _MData_ vIvL_f    = {       ivL_f,       ivL_f,       ivL_f,       ivL_f,       ivL_f,       ivL_f,       ivL_f,       ivL_f,
				      ivL_f,       ivL_f,       ivL_f,       ivL_f,       ivL_f,       ivL_f,       ivL_f,       ivL_f };
#elif   defined(__AVX__)
constexpr _MData_ oPif      = {     1./M_PI,     1./M_PI,     1./M_PI,     1./M_PI,     1./M_PI,     1./M_PI,     1./M_PI,     1./M_PI };
constexpr _MData_ zeroNegf  = {       -0.0f,       -0.0f,       -0.0f,       -0.0f,       -0.0f,       -0.0f,       -0.0f,       -0.0f };
constexpr _MData_ fHlf       = {        0.5f,        0.5f,        0.5f,        0.5f,        0.5f,        0.5f,        0.5f,        0.5f };
#ifdef	__FMA__
constexpr _MData_ TriMaxf   = {         1e7,         1e7,         1e7,         1e7,         1e7,         1e7,         1e7,         1e7 };
#else
constexpr _MData_ TriMaxf   = {         1e5,         1e5,         1e5,         1e5,         1e5,         1e5,         1e5,         1e5 };
#endif
constexpr _MData_ fInf      = {       Inf_f,       Inf_f,       Inf_f,       Inf_f,       Inf_f,       Inf_f,       Inf_f,       Inf_f };
constexpr _MData_ fNan      = {       Nan_f,       Nan_f,       Nan_f,       Nan_f,       Nan_f,       Nan_f,       Nan_f,       Nan_f };
constexpr _MData_ PiAf      = {       PiA_f,       PiA_f,       PiA_f,       PiA_f,       PiA_f,       PiA_f,       PiA_f,       PiA_f };
constexpr _MData_ PiBf      = {       PiB_f,       PiB_f,       PiB_f,       PiB_f,       PiB_f,       PiB_f,       PiB_f,       PiB_f };
constexpr _MData_ PiCf      = {       PiC_f,       PiC_f,       PiC_f,       PiC_f,       PiC_f,       PiC_f,       PiC_f,       PiC_f };
constexpr _MData_ PiDf      = {       PiD_f,       PiD_f,       PiD_f,       PiD_f,       PiD_f,       PiD_f,       PiD_f,       PiD_f };
constexpr _MData_ hPiAf     = {  0.5f*PiA_f,  0.5f*PiA_f,  0.5f*PiA_f,  0.5f*PiA_f,  0.5f*PiA_f,  0.5f*PiA_f,  0.5f*PiA_f,  0.5f*PiA_f };
constexpr _MData_ hPiBf     = {  0.5f*PiB_f,  0.5f*PiB_f,  0.5f*PiB_f,  0.5f*PiB_f,  0.5f*PiB_f,  0.5f*PiB_f,  0.5f*PiB_f,  0.5f*PiB_f };
constexpr _MData_ hPiCf     = {  0.5f*PiC_f,  0.5f*PiC_f,  0.5f*PiC_f,  0.5f*PiC_f,  0.5f*PiC_f,  0.5f*PiC_f,  0.5f*PiC_f,  0.5f*PiC_f };
constexpr _MData_ hPiDf     = {  0.5f*PiD_f,  0.5f*PiD_f,  0.5f*PiD_f,  0.5f*PiD_f,  0.5f*PiD_f,  0.5f*PiD_f,  0.5f*PiD_f,  0.5f*PiD_f };
constexpr _MData_ s0f       = {        s0_f,        s0_f,        s0_f,        s0_f,        s0_f,        s0_f,        s0_f,        s0_f };
constexpr _MData_ s1f       = {        s1_f,        s1_f,        s1_f,        s1_f,        s1_f,        s1_f,        s1_f,        s1_f };
constexpr _MData_ s2f       = {        s2_f,        s2_f,        s2_f,        s2_f,        s2_f,        s2_f,        s2_f,        s2_f };
constexpr _MData_ s3f       = {        s3_f,        s3_f,        s3_f,        s3_f,        s3_f,        s3_f,        s3_f,        s3_f };

constexpr _MData_ vTh1f     = {       th1_f,       th1_f,       th1_f,       th1_f,       th1_f,       th1_f,       th1_f,       th1_f };
constexpr _MData_ vTh2f     = {       th2_f,       th2_f,       th2_f,       th2_f,       th2_f,       th2_f,       th2_f,       th2_f };
constexpr _MData_ vL1f      = {        L1_f,        L1_f,        L1_f,        L1_f,        L1_f,        L1_f,        L1_f,        L1_f };
constexpr _MData_ vL2f      = {       -L2_f,       -L2_f,       -L2_f,       -L2_f,       -L2_f,       -L2_f,       -L2_f,       -L2_f };
constexpr _MData_ vA1f      = {        A1_f,        A1_f,        A1_f,        A1_f,        A1_f,        A1_f,        A1_f,        A1_f };
constexpr _MData_ vA2f      = {        A2_f,        A2_f,        A2_f,        A2_f,        A2_f,        A2_f,        A2_f,        A2_f };
constexpr _MData_ vIvL_f    = {       ivL_f,       ivL_f,       ivL_f,       ivL_f,       ivL_f,       ivL_f,       ivL_f,       ivL_f };
#else
constexpr _MData_ oPif      = {     1./M_PI,     1./M_PI,     1./M_PI,     1./M_PI };
constexpr _MData_ zeroNegf  = {       -0.0f,       -0.0f,       -0.0f,       -0.0f };
constexpr _MData_ fHlf       = {        0.5f,        0.5f,        0.5f,        0.5f };
constexpr _MData_ TriMaxf   = {         1e5,         1e5,         1e5,         1e5 };
constexpr _MData_ fInf      = {       Inf_f,       Inf_f,       Inf_f,       Inf_f };
constexpr _MData_ PiAf      = {       PiA_f,       PiA_f,       PiA_f,       PiA_f };
constexpr _MData_ PiBf      = {       PiB_f,       PiB_f,       PiB_f,       PiB_f };
constexpr _MData_ PiCf      = {       PiC_f,       PiC_f,       PiC_f,       PiC_f };
constexpr _MData_ PiDf      = {       PiD_f,       PiD_f,       PiD_f,       PiD_f };
constexpr _MData_ hPiAf     = {  0.5f*PiA_f,  0.5f*PiA_f,  0.5f*PiA_f,  0.5f*PiA_f };
constexpr _MData_ hPiBf     = {  0.5f*PiB_f,  0.5f*PiB_f,  0.5f*PiB_f,  0.5f*PiB_f };
constexpr _MData_ hPiCf     = {  0.5f*PiC_f,  0.5f*PiC_f,  0.5f*PiC_f,  0.5f*PiC_f };
constexpr _MData_ hPiDf     = {  0.5f*PiD_f,  0.5f*PiD_f,  0.5f*PiD_f,  0.5f*PiD_f };
constexpr _MData_ s0f       = {        s0_f,        s0_f,        s0_f,        s0_f };
constexpr _MData_ s1f       = {        s1_f,        s1_f,        s1_f,        s1_f };
constexpr _MData_ s2f       = {        s2_f,        s2_f,        s2_f,        s2_f };
constexpr _MData_ s3f       = {        s3_f,        s3_f,        s3_f,        s3_f };

constexpr _MData_ vTh1f     = {       th1_f,       th1_f,       th1_f,       th1_f };
constexpr _MData_ vTh2f     = {       th2_f,       th2_f,       th2_f,       th2_f };
constexpr _MData_ vL1f      = {        L1_f,        L1_f,        L1_f,        L1_f };
constexpr _MData_ vL2f      = {       -L2_f,       -L2_f,       -L2_f,       -L2_f };
constexpr _MData_ vA1f      = {        A1_f,        A1_f,        A1_f,        A1_f };
constexpr _MData_ vA2f      = {        A2_f,        A2_f,        A2_f,        A2_f };
constexpr _MData_ vIvL_f    = {       ivL_f,       ivL_f,       ivL_f,       ivL_f };
#endif

constexpr std::array<float, 32> sTrail_f = { St00_f, St01_f, St02_f, St03_f, St04_f, St05_f, St06_f, St07_f,
       					     St08_f, St09_f, St10_f, St11_f, St12_f, St13_f, St14_f, St15_f,
       					     St16_f, St17_f, St18_f, St19_f, St20_f, St21_f, St22_f, St23_f,
       					     St24_f, St25_f, St26_f, St27_f, St28_f, St29_f, St30_f, St31_f };

constexpr std::array<float, 32> sLead_f  = { Sl00_f, Sl01_f, Sl02_f, Sl03_f, Sl04_f, Sl05_f, Sl06_f, Sl07_f,
       					     Sl08_f, Sl09_f, Sl10_f, Sl11_f, Sl12_f, Sl13_f, Sl14_f, Sl15_f,
       					     Sl16_f, Sl17_f, Sl18_f, Sl19_f, Sl20_f, Sl21_f, Sl22_f, Sl23_f,
       					     Sl24_f, Sl25_f, Sl26_f, Sl27_f, Sl28_f, Sl29_f, Sl30_f, Sl31_f };


/*	Sleef	*/
inline void printiVar(_MInt_ d, const char *name) {

	printf ("%s", name);
#if	defined(__AVX512F__)
	int r[16] __attribute((aligned(64)));
	opCode(store_si512, r, d);
	for (int i=0; i<16; i++)
#elif	defined(__AVX__)
	int r[8] __attribute((aligned(32)));
	opCode(store_si256, ((_MInt_ *)r), d);
	for (int i=0; i<8; i++)
#else
	int r[4] __attribute((aligned(16)));
	opCode(store_si128, ((_MInt_ *)r), d);
	for (int i=0; i<4; i++)
#endif
		printf(" %d", r[i]);
	printf("\n");
}

inline void printsVar(_MData_ d, const char *name) {

	printf ("%s", name);
#if	defined(__AVX512F__)
	for (int i=0; i<16; i++)
#elif	defined(__AVX__)
	for (int i=0; i<8; i++)
#else
	for (int i=0; i<4; i++)
#endif
		printf(" %f", d[i]);
	printf("\n");
}

#endif

#undef	_MData_
#undef	_PREFIX
#undef opCode_P
#undef opCode_N
#undef opCode

#endif
