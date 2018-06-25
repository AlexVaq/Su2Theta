#ifndef	__RANDOMGEN
	#define	__RANDOMGEN
	#include <vector>
	#include <random>
	#include <memory>

	namespace	Su2Rand {

		void	initRandom ();

		template<typename Float>
		class	RandomGen {
			private:

			std::vector<int>			sd;
			std::vector<std::mt19937_64>		mt64;
			std::uniform_real_distribution<Float>	uni;

			public:

				RandomGen();
				RandomGen(int  seed);

			Float	rand();
		};

		extern std::shared_ptr<RandomGen<double>> myRNG;

		inline auto	genRand () { return (*myRNG).rand(); }
	}
#endif
