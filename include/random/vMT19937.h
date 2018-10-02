#ifndef	__VMT19937
	#define	__VMT19937
	#include <random>
	#include <memory>

	#include "enumFields.h"
	#include "comms/comms.h"
	#include "utils/utils.h"

	#include "simd/simd.h"

	#include <omp.h>

	using namespace	Su2Enum;

	namespace	Su2Rand {

		constexpr int    mtSl32		= 18;
		constexpr int    mtSr32		= 11;
		constexpr int    mtPos32	= 122;
		constexpr int    mtSize32	= 156;
		constexpr int    mtPos64	= 117;
		constexpr int    mtSize64	= 191;
		constexpr int    mtSl64		= 19;
		constexpr int    mtSr64		= 12;
		constexpr uint   mtInit		= 1812433253UL;
		constexpr uint64 mtFix1		= 0x90014964b32f4329;
		constexpr uint64 mtFix2		= 0x3b8d12ac548a7c7a;
		constexpr uint64 mtPcd1		= 0x3d84e1ac0dc82880;
		constexpr uint64 mtPcd2		= 0x0000000000000001;
		constexpr uint	 mtPcf1		= 0x00000001;
		constexpr uint	 mtPcf2		= 0x13c9e684;
		constexpr uint   mtLowMask32	= 0x007fffff;
		constexpr uint64 mtLowMask64	= 0x000fffffffffffff;
		constexpr uint   mtHighConst32	= 0x3f800000;
		constexpr uint64 mtHighConst64	= 0x3ff0000000000000;
		const	  Simd::Simd_f mtMsk32(0xdfffffefU, 0xddfecb7fU, 0xbffaffffU, 0xbffffff6U);
		const	  Simd::Simd_d mtMsk64(0x000ffafffffffb3fU,      0x000ffdfffc90fffdU);

		const	  Simd::Simd_f mtLow32 (mtLowMask32,   mtLowMask32,   mtLowMask32,   mtLowMask32);
		const	  Simd::Simd_f mtHigh32(mtHighConst32, mtHighConst32, mtHighConst32, mtHighConst32);

		void	vInitRandom ();

		template<class T>
		class	RNG {
			private:

			//T	*state;	// Not a big deal for Simd_f and this way we have a single function
			T	state[mtSize64+1]  __attribute__ ((aligned (Simd::sAlign)));
			int	idx;

			void	initRNG (uint seed);
			void	genNext ();

			public:

				RNG ()		{
//				alignAlloc((void **) &state, Simd::sAlign, (mtSize64+1)*sizeof(T));
				initRNG(0);
			}
				RNG (uint seed)	{
//				alignAlloc((void **) &state, Simd::sAlign, (mtSize64+1)*sizeof(T));
				initRNG(seed);
			}

//				~RNG () {
//				trackFree(state);
//			}

			void	Seed(uint seed)	{ initRNG(seed); }

			T	operator()();
		};

		template<class T>
		class	vRandomGen {
			private:

			std::vector<uint>			sd;
			std::vector<RNG<T>>			mt;

			public:

				vRandomGen() {

				sd.resize(Su2Comms::nThreads);
				mt.resize(Su2Comms::nThreads);

				std::random_device seed;			// Totally random seed coming from memory garbage

				for (int i=0; i<Su2Comms::nThreads; i++)
					sd[i] = i;//seed()*(1 + Su2Comms::rank);

				#pragma omp parallel default(shared)
				{
					int nThread = omp_get_thread_num();
					mt[nThread].Seed(sd[nThread]);	// Vectorized Mersenne-Twister, independent per thread
				}
			}

				vRandomGen(int seed) {

				sd.resize(Su2Comms::nThreads);
				mt.resize(Su2Comms::nThreads);

				for (int i=0; i<Su2Comms::nThreads; i++)
					sd[i] = seed*(1 + Su2Comms::rank);

				#pragma omp parallel default(shared)
				{
					int nThread = omp_get_thread_num();
					mt[nThread].Seed(sd[nThread]);	// Vectorized Mersenne-Twister, independent per thread
				}
			}

			T	operator()() {
				int nThread = omp_get_thread_num();
				return	mt[nThread]();
			}
		};

		template<class T>
		T	vGenRand();
	}


#endif
