#ifndef __ENUM_FIELDS
	#define __ENUM_FIELDS
	#include <cstdint>
	#include <string>
	#include <mpi.h>

	typedef	unsigned int	uint;
	typedef	uint64_t	uint64;

	namespace	Su2Enum {

		#if   defined __AVX512F__
		const std::string SystemVec{"Avx512"};
		#elif defined __AVX2__
		const std::string SystemVec{"Avx2"};
		#elif defined __AVX__
		const std::string SystemVec{"Avx"};
		#else
		const std::string SystemVec{"Sse4"};
		#endif

		typedef	enum	Device_s {
			DeviceCpu,
			DeviceGpu,
			DeviceGpuManaged,
		}	Device;

		typedef	enum	Precision_s {
			PrecDouble = 8,
			PrecSingle = 4,
		}	Precision;

		typedef enum    LogLevel_s
		{
			MsgLog   = 1048576,
			DebugLog = 2097152,
			ErrorLog = 4194304,
		}	LogLevel;

		typedef	enum	ColorKind_s
		{
			EvenOdd	= 2,
			Colored	= 32,
		}	ColorKind;

		typedef enum    ParityType_s
		{
			ParityEven = 0,
			ParityOdd  = 1,
		}	ParityType;

		typedef enum    Parity2Type_s
		{
			Parity2Even = 0,
			Parity2Odd  = 1,
		}	Parity2Type;

		typedef enum    ParityColor_s
		{
			BlackEven     = 0,
			WhiteEven     = 1,
			DarkRedEven   = 2,
			RedEven       = 3,
			CoralEven     = 4,
			OrangeEven    = 5,
			SandEven      = 6,
			YellowEven    = 7,
			GreenEven     = 8,
			LimeEven      = 9,
			TealEven      = 10,
			SteelEven     = 11,
			BlueEven      = 12,
			SlateEven     = 13,
			VioletEven    = 14,
			OrchidEven    = 15,
			BlackOdd      = 16,
			WhiteOdd      = 17,
			DarkRedOdd    = 18,
			RedOdd        = 19,
			CoralOdd      = 20,
			OrangeOdd     = 21,
			SandOdd       = 22,
			YellowOdd     = 23,
			GreenOdd      = 24,
			LimeOdd       = 25,
			TealOdd       = 26,
			SteelOdd      = 27,
			BlueOdd       = 28,
			SlateOdd      = 29,
			VioletOdd     = 30,
			OrchidOdd     = 31,
		}	ParityColor;

		typedef enum    CommOperation_s
		{
			CommSend,
			CommRecv,
			CommSdRv,
			CommWait,
		}	CommOperation;

		typedef enum    LogMpi_s
		{
			AllRanks,
			ZeroRank,
		}	LogMpi;

		typedef enum    VerbosityLevel_s
		{
			VerbSilent = 0,
			VerbNormal = 1,
			VerbHigh   = 2,
		}	VerbosityLevel;

		typedef enum    FindType_s {
			FindMax,
			FindMin,
		}	FindType;

		typedef enum    AllocType_s
		{
			AllocTrack = 0,
			AllocAlign = 1,
		}	AllocType;

		typedef enum    ProfType_s
		{
			ProfGen,
			ProfPlaq,
			ProfAction,
			ProfHB,
			ProfOvR,
			ProfMetro,
			ProfQCharge,
			ProfTuner,
			ProfFold,
			ProfHdf5,
		}	ProfType;

		typedef	enum	Permutation_s
		{
			NoPermutation   = 0b0000,
			PermutationX    = 0b0001,
			PermutationY    = 0b0010,
			PermutationZ    = 0b0100,
			PermutationT    = 0b1000,
			PermutationXY   = 0b0011,
			PermutationXZ   = 0b0101,
			PermutationXT   = 0b1001,
			PermutationYZ   = 0b0110,
			PermutationYT   = 0b1010,
			PermutationZT   = 0b1100,
			PermutationXYZ  = 0b0111,
			PermutationXYT  = 0b1011,
			PermutationXZT  = 0b1101,
			PermutationYZT  = 0b1110,
			PermutationXYZT = 0b1111,
		}	Permutation;

#ifdef	__NVCC__
	#define	Attr	inline constexpr __host__ __device__
#else
	#define	Attr	inline constexpr
#endif
		template<typename enumFlag>
		Attr enumFlag  operator &  (enumFlag  lhs, const enumFlag rhs) { return static_cast<enumFlag>(static_cast<int>(lhs) & static_cast<const int>(rhs)); }
		template<typename enumFlag>
		Attr enumFlag& operator &= (enumFlag &lhs, const enumFlag rhs) { lhs  = static_cast<enumFlag>(static_cast<int>(lhs) & static_cast<const int>(rhs)); return lhs; }
		template<typename enumFlag>
		Attr enumFlag  operator |  (enumFlag  lhs, const enumFlag rhs) { return static_cast<enumFlag>(static_cast<int>(lhs) | static_cast<const int>(rhs)); }
		template<typename enumFlag>
		Attr enumFlag& operator |= (enumFlag &lhs, const enumFlag rhs) { lhs  = static_cast<enumFlag>(static_cast<int>(lhs) | static_cast<const int>(rhs)); return lhs; }
		template<typename enumFlag>
		Attr enumFlag  operator ^  (enumFlag  lhs, const enumFlag rhs) { return static_cast<enumFlag>(static_cast<int>(lhs) ^ static_cast<const int>(rhs)); }
		template<typename enumFlag>
		Attr enumFlag& operator ^= (enumFlag &lhs, const enumFlag rhs) { lhs  = static_cast<enumFlag>(static_cast<int>(lhs) ^ static_cast<const int>(rhs)); return lhs; }
#undef	Attr
	}
#endif
