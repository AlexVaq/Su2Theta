#ifndef __ENUM_FIELDS
	#define __ENUM_FIELDS
	#include <mpi.h>

	namespace	Su2Enum {

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

		typedef enum    ParityType_s
		{
			ParityEven = 0,
			ParityOdd  = 1,
		}	ParityType;

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
			ProfQCharge,
			ProfTuner,
			ProfPlaq,
			ProfFold,
			ProfHdf5,
		}	ProfType;
	}
#endif
