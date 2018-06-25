#include <memory>
#include "utils/logger.h"

using namespace Su2Enum;

namespace Su2Log {
	std::shared_ptr<Logger> myLog;
	const char	levelTable[3][16] = { " Msg ", "Debug", "Error" };
}

void	createLogger(const int index, const LogMpi logMpi, const VerbosityLevel verbosity) {
	static bool	initLog = false;

	if	(initLog == false) {
		Su2Log::myLog = std::make_shared<Su2Log::Logger>(index, logMpi, verbosity);
	}	else	{
		LogMsg(VerbHigh, "Logger already initialized");
	}
}


