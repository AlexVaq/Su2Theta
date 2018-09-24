#ifndef	_TUNABLE_
	#define	_TUNABLE_
	#include <string>
	#include <utils/profiler.h>

	using namespace Su2Prof;

	class	Tunable {
		protected:

		std::string	name;

		double		gFlops;
		double		gBytes;

		unsigned int	xBlock;
		unsigned int	yBlock;
		unsigned int	zBlock;
		unsigned int	tBlock;

		unsigned int	xBest;
		unsigned int	yBest;
		unsigned int	zBest;
		unsigned int	tBest;

		unsigned int	xMax;
		unsigned int	yMax;
		unsigned int	zMax;
		unsigned int	tMax;

		bool		isTuned;
		bool		isGpu;

		public:

				Tunable() noexcept : gFlops(0.), gBytes(0.), xBlock(0), yBlock(0), zBlock(0), tBlock(0), xBest(0), yBest(0), zBest(0), tBest(0),
						     isTuned(false), isGpu(false), name("") {}

		inline	double		GFlops () const noexcept { return gFlops; }
		inline	double		GBytes () const noexcept { return gBytes; }

		inline	unsigned int	BlockX () const noexcept { return xBlock; }
		inline	unsigned int	BlockY () const noexcept { return yBlock; }
		inline	unsigned int	BlockZ () const noexcept { return zBlock; }
		inline	unsigned int	BlockT () const noexcept { return tBlock; }

		inline	bool		IsTuned() const noexcept { return isTuned;  }
			void		UnTune ()       noexcept { isTuned = false; }
			void		Tune   ()       noexcept { isTuned = true;  }

		inline	unsigned int	TunedBlockX () const noexcept { return xBest; }
		inline	unsigned int	TunedBlockY () const noexcept { return yBest; }
		inline	unsigned int	TunedBlockZ () const noexcept { return zBest; }
		inline	unsigned int	TunedBlockT () const noexcept { return tBest; }

		inline	unsigned int	MaxBlockX () const noexcept { return xMax; }
		inline	unsigned int	MaxBlockY () const noexcept { return yMax; }
		inline	unsigned int	MaxBlockZ () const noexcept { return zMax; }
		inline	unsigned int	MaxBlockT () const noexcept { return tMax; }

		inline	size_t		TotalThreads() const noexcept { return xBlock*yBlock*zBlock*tBlock; }

		inline	unsigned int	SetBlockX (unsigned int bSize) noexcept { xBlock = bSize; }
		inline	unsigned int	SetBlockY (unsigned int bSize) noexcept { yBlock = bSize; }
		inline	unsigned int	SetBlockZ (unsigned int bSize) noexcept { zBlock = bSize; }
		inline	unsigned int	SetBlockT (unsigned int bSize) noexcept { tBlock = bSize; }

		inline	void		UpdateBestBlock() noexcept { xBest  = xBlock; yBest  = yBlock; zBest  = zBlock; tBest  = tBlock; }
		inline	void		SetBestBlock()    noexcept { xBlock = xBest;  yBlock = yBest;  zBlock = zBest;  tBlock = tBest;  }

		void		AdvanceBlockSize() noexcept {

			if (isGpu) {
/*				do {
					if (xBlock < xMax) {
						do {
							xBlock++;
						}	while ((xSize % xBlock) != 0);
					} else {
						xBlock = 8;

						if (yBlock < yMax) {
							do {
								yBlock++;
							}	while ((ySize % yBlock) != 0);
						} else {
							isTuned = true;
						}
					}
				}	while (!isTuned && TotalThreads() > maxThreadsPerBlock());
				*/
			} else {
				if (xBlock < xMax) {
					do {
						xBlock+=2;
					}	while ((xMax % xBlock) != 0);
				} else {
					xBlock = 4;

					if (yBlock < yMax) {
						do {
							yBlock+=2;
						}	while ((yMax % yBlock) != 0);
					} else {
						yBlock = 2;

						if (zBlock < zMax) {
							do {
								zBlock+=2;
							}	while ((zMax % zBlock) != 0);
						} else {
							zBlock = 2;

							if (tBlock < tMax) {
								do {
									tBlock+=2;
								}	while ((tMax % tBlock) != 0);
							} else {
								isTuned = true;
							}
						}
					}
				}
			}
		}	

		std::string	Name   () const noexcept { return name; }

		void		reset  ()                     { gFlops = 0.; gBytes = 0.; }
		void		add    (double GF, double GB) { gFlops += GF; gBytes += GB; }

		void		SetName   (std::string  newName) { name = newName; }
		void		SetName   (const char * newName) { name.assign(newName); }
		void		AppendName(std::string  appName) { name += appName; }
		void		AppendName(const char * appName) { name += std::string(appName); }

		void		ResetBlockSize(bool gpu = false) {
			if (!isGpu) {
				xBest = xBlock = 4;
				yBest = yBlock = 2;
				zBest = zBlock = 2;
				tBest = tBlock = 2;
			} else {
				xBest = xBlock = 8;
				yBest = yBlock = 1;
				zBest = zBlock = 1;
			}
		}

		void		InitBlockSize(size_t vL[4], bool gpu = false) {
			size_t lV = SIZE_MAX;

			isGpu = gpu;

			if (!isGpu) {
				xBest = xBlock = xMax = vL[0];
				yBest = yBlock = yMax = vL[1];
				zBest = zBlock = zMax = vL[2];
				tBest = tBlock = tMax = vL[3];
//			} else {
//				auto xTmp = maxThreadsPerDim(0); 
//				auto yTmp = maxThreadsPerDim(1);
//				auto zTmp = 1;
//
//				lV = maxThreadsPerBlock();
//
//				xMax = (Lx*Lx > xTmp) ? xTmp : Lx*Lx;
//				yMax = (Lz > yTmp) ? yTmp : Lz;
//				zMax = 1;
//
//				xSize = Lx*Lx;
//				ySize = Lz;
//				zSize = 1;
//
//				if (yTmp*maxGridSize(2) < Lz)
//					LogError("Error: not enough threads on gpu to accomodate z-dimension");
//
//				xBest = xBlock = 8;
//				yBest = yBlock = 1;
//				zBest = zBlock = 1;
			}

			isTuned = false;
		}

		void		*backup;

		virtual	void	SaveState	() {};
		virtual	void	RestoreState	() {};
		virtual	void	Run		() {};
		virtual	void	Reset		() {};
	};

	namespace	Su2Tune {
		void	Tune	(Tunable &func);
	}
	
#endif
