#ifndef	__LOGGER__
	#define	__LOGGER__

	#include <string>
	#include <cstring>
	#include <sstream>
	#include <fstream>
	#include <iomanip>
	#include <vector>
	#include <chrono>
	#include <memory>
	#include <cstddef>
	#include <algorithm>
	#include <thread>

	#include <cstdarg>
	#include <cstring>
	#include <sys/stat.h>
	#include <unistd.h>

	#include <omp.h>
	#include <mpi.h>

	#include "enumFields.h"
	#include "comms/comms.h"

	using namespace	Su2Enum;

	namespace Su2Log {

		constexpr long long int	logFreq	= 50000;
		constexpr int		MsgSize	= 2048;
		constexpr size_t 	basePack = sizeof(ptrdiff_t)*5;

		extern	const char	levelTable[3][16];

		class	Msg {
			private:
				std::chrono::time_point<std::chrono::high_resolution_clock> timestamp;	// Timestamp

				int		tIdx;							// OMP thread
				int		mRnk;							// MPI rank
				int		size;							// Message size
				std::string	data;							// Log message
				LogLevel	logLevel;						// Level of logging (info, debug or error)

				char		packed[MsgSize];					// For serialization and MPI

				mutable MPI_Request req;

			public:
					Msg(LogLevel logLevel, const int tIdx, const char * format, ...) noexcept : logLevel(logLevel), tIdx(tIdx) {
					char buffer[MsgSize - basePack];
					va_list args;
					va_start (args, format);
					size = vsnprintf (buffer, MsgSize - 1 - basePack, format, args);
					va_end (args);

					data.assign(buffer, size);

					timestamp = std::chrono::high_resolution_clock::now();
					mRnk = Su2Comms::rank;
				}

					 Msg(const Msg &myMsg) noexcept : tIdx(myMsg.tIdx), size(myMsg.size), data(myMsg.data), logLevel(myMsg.logLevel),
									  mRnk(myMsg.mRnk), timestamp(myMsg.timestamp), req(MPI_REQUEST_NULL) {};
					 Msg(Msg &&myMsg)      noexcept : tIdx(std::move(myMsg.tIdx)), size(std::move(myMsg.size)), data(std::move(myMsg.data)), mRnk(std::move(myMsg.mRnk)),
									  logLevel(std::move(myMsg.logLevel)), timestamp(myMsg.timestamp), req(MPI_REQUEST_NULL) {};
					 Msg(void *vPack)      noexcept : tIdx(static_cast<ptrdiff_t*>(vPack)[1]), size(static_cast<ptrdiff_t*>(vPack)[3]),
									  logLevel((LogLevel) static_cast<ptrdiff_t*>(vPack)[2]), mRnk(static_cast<ptrdiff_t*>(vPack)[4]) {
					auto mTime = static_cast<ptrdiff_t*>(vPack)[0];
					timestamp  = std::chrono::time_point<std::chrono::high_resolution_clock>(std::chrono::microseconds(mTime));

					req = MPI_REQUEST_NULL;
					char *msgData = static_cast<char*>(vPack) + basePack;
					data.assign(msgData, size);
				};

					~Msg() noexcept {};

				Msg&	operator=(const Su2Log::Msg& msg) {
					tIdx = msg.tIdx;
					size = msg.size;
					data = msg.data;
					mRnk = msg.mRnk;

					logLevel  = msg.logLevel;
					timestamp = msg.timestamp;

					return	(*this);
				}

				inline	long long int	time(std::chrono::time_point<std::chrono::high_resolution_clock> start) const {
					return std::chrono::duration_cast<std::chrono::microseconds> (timestamp - start).count();
				}

				inline	int		thread() const { return tIdx; }
				inline	int		rank()   const { return mRnk; }
				inline	std::string	msg()    const { return data; }
				inline	const LogLevel	level()  const { return logLevel; }

				inline	char*		pack() {
					ptrdiff_t *sPack = static_cast<ptrdiff_t*>(static_cast<void*>(packed));
					sPack[0] = std::chrono::time_point_cast<std::chrono::microseconds> (timestamp).time_since_epoch().count();
					sPack[1] = tIdx;
					sPack[2] = (ptrdiff_t) logLevel;
					sPack[3] = data.length();
					sPack[4] = mRnk;

					memcpy (packed+basePack, data.data(), data.length());

					return	packed;
				}

				inline	MPI_Request&	mpiReq() { return req; }
		};

		class	Logger {
			private:
				std::chrono::time_point<std::chrono::high_resolution_clock> logStart;
				std::ofstream		oFile;
				const LogMpi		mpiType;
				std::vector<Msg>	msgStack;
				const VerbosityLevel	verbose;
				bool			logRunning;
				bool			logWriting;

				void	printMsg	(const Msg &myMsg) noexcept {
					oFile << std::setw(11) << myMsg.time(logStart)/1000 << "ms: Logger level[" << std::right << std::setw(5) << levelTable[myMsg.level()>>21] << "]"
					      << " Rank " << std::setw(4)  << myMsg.rank()+1 << "/" << Su2Comms::size << " - Thread " << std::setw(3) << myMsg.thread()+1 << "/"
					      << omp_get_num_threads() << " ==> " << myMsg.msg() << std::endl;
				}

				/* We only allow thread 0 to write to disk, but any other thread can put messages in the stack		*/
				/* The stack is flushed if there is an error on any thread because the variable mustFlush is shared	*/
				void	flushMsg	() noexcept {
					if (omp_get_thread_num() != 0 || logWriting == true)
						return;

					logWriting = true;
					auto it = msgStack.cbegin();

					if (it != msgStack.cend()) {
						printMsg(*it);
						msgStack.erase(it);
					}
					logWriting = false;
				}

				void	flushStack	() noexcept {
					if (omp_get_thread_num() != 0)
						return;

					std::sort(msgStack.begin(), msgStack.end(), [logStart = logStart](Msg a, Msg b) { return (a.time(logStart) < b.time(logStart)); } );

					for (auto it = msgStack.cbegin(); it != msgStack.cend(); it++)
						printMsg(*it);

					msgStack.clear();
				}

				void	flushDisk	() {
					if (omp_get_thread_num() != 0)
						return;
					oFile.flush();
				}

				bool	getMpiMsg	(const LogLevel level) {
					if (omp_get_thread_num() != 0)
						return	false;

					bool msgPending = false;

					int flag = 0;
					MPI_Status status;

					// Get all the messages of a particular log level
					do {
						MPI_Iprobe(MPI_ANY_SOURCE, level, MPI_COMM_WORLD, &flag, &status);

						if (flag) {
							char packed[MsgSize];
							int  mSize;
							auto srcRank = status.MPI_SOURCE;

							// Get message
							MPI_Get_count(&status, MPI_CHAR, &mSize);
							MPI_Recv(packed, mSize, MPI_CHAR, srcRank, level, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							// Put the message in the stack
							while (logWriting == true) {}
							msgStack.emplace_back(static_cast<void*>(packed));
							msgPending = true;
						}
					}	while (flag);

					return	msgPending;
				}

			public:
				 Logger(const int index, const LogMpi mpiType, const VerbosityLevel verbosity) : mpiType(mpiType), verbose(verbosity) {

					bool			test;
					struct stat		buffer;
					std::string		base("Su2.log.");
					std::stringstream	ss;

					int			idx = index - 1;

					// Let's force a sync before recording the starting time, so all the ranks more or less agree on this
					Su2Comms::commSync();
					logStart = std::chrono::high_resolution_clock::now();

					if (Su2Comms::rank == 0) {
						logRunning = true;
						do {
							idx++;
							ss.str("");
							ss << base << idx;
							test = (stat (ss.str().c_str(), &buffer) == 0) ? true : false;
						}	while (test);

						oFile.open(ss.str().c_str(), std::ofstream::out);
						banner();

						std::thread([&]() {
							while (logRunning) {
								std::this_thread::sleep_for(std::chrono::microseconds(logFreq));
								logWriting = true;
								flushLog();
								logWriting = false;
							}
						}).detach();
					}
					Su2Comms::commSync();

				}

				// Stops logger in preparation for an MPI shutdown
				void	stop		() { logRunning = false; usleep(logFreq); }

				// Receives pending MPI messages and flushes them to disk
				void	flushLog	() {
					if (omp_get_thread_num() != 0)
						return;

					// Get all the messages
					if (mpiType == AllRanks) {
						getMpiMsg (MsgLog);
						getMpiMsg (ErrorLog);
					}

					flushStack();
					flushDisk();
				}

				~Logger() { int noMpi; MPI_Finalized(&noMpi); if (noMpi == 0) flushLog(); if (Su2Comms::rank==0) { oFile.close(); } }

				auto	runTime() {
					auto	cTime = std::chrono::high_resolution_clock::now();
					auto	dTime = std::chrono::duration_cast<std::chrono::microseconds> (cTime - logStart).count();

					return	dTime;
				}

				template<typename... Fargs>
				void	operator()(LogLevel level, const char * file, const int line, const char * format, Fargs... vars)
				{
					static bool       mustFlush = false;
					static const int  myRank    = Su2Comms::rank;
					static const int  nSplit    = Su2Comms::size;
					static const bool mpiLogger = ((nSplit > 1) && (mpiType == AllRanks));

					if (mpiType == ZeroRank && myRank != 0)
						return;

					switch	(level) {

						case	MsgLog:
						{
							// We push the messages in the stack and we flush them later
							msgStack.push_back(std::move(Msg(level, omp_get_thread_num(), format, vars...)));

//							if	(std::chrono::duration_cast<std::chrono::microseconds> (std::chrono::high_resolution_clock::now() - logStart).count() > logFreq)
//								mustFlush = true;
							break;
						}

						case	DebugLog:
						{
							Msg	msg(level, omp_get_thread_num(), format, vars...);
							printMsg(msg);										// We directly write to disk
							flushDisk();
							break;
						}

						case	ErrorLog:
						{
							auto thisThread = omp_get_thread_num();
							// We add a message to the stack and immediately flush the whole message stack
							msgStack.push_back(std::move(Msg(level, thisThread, "Error in file %s line %d", file, line)));
							msgStack.push_back(std::move(Msg(level, thisThread, format, vars...)));
							mustFlush = true;
							break;
						}
					}

					if (mpiLogger) {	// If we are using MPI
						if (myRank == 0) {
							// Get all the pending error messages. If found any, prepare to flush the stack
							if (getMpiMsg (ErrorLog))	// We use this if not to overwrite mustFlush
								mustFlush = true;
						} else {
							/* If loglevel is ERROR, we block comms until all the pending messages are sent and  we clear the stack */
							if (level == ErrorLog) {
								for (auto it = msgStack.begin(); it != msgStack.end();) {
									MPI_Request& req = it->mpiReq();
									if (req == MPI_REQUEST_NULL)
										MPI_Isend(it->pack(), it->msg().length() + basePack, MPI_CHAR, 0, level, MPI_COMM_WORLD, &req);
									MPI_Wait(&req, MPI_STATUS_IGNORE);
									it = msgStack.erase(it);
								}
								msgStack.clear();
							} else {
								/* As the messages are sent, we remove them from the stack */
								for (auto it = msgStack.begin(); it != msgStack.end();) {
									int flag = 0;
									MPI_Request& req = it->mpiReq();
									if (req == MPI_REQUEST_NULL) {
										MPI_Isend(it->pack(), it->msg().length() + basePack, MPI_CHAR, 0, level, MPI_COMM_WORLD, &req);
										it++;
									} else {
										MPI_Test(&req, &flag, MPI_STATUS_IGNORE);

										if (flag)
											it = msgStack.erase(it);
										else
											it++;
									}
								}
							}

						}
					}

					if (mustFlush && myRank == 0) {
						if (mpiLogger) {	// If we are using MPI
							// Get all the standard messages
							getMpiMsg(MsgLog);
						}

						flushStack();
						flushDisk();
						mustFlush = false;
					}
				}

				void	banner		() {
					(*this)(MsgLog, nullptr, 0, "Su2 logger started");
				}

				const VerbosityLevel	Verbosity	() const { return	verbose; }
		};

		extern std::shared_ptr<Logger> myLog;
	};

	void	createLogger(const int index, const LogMpi logMpi, const VerbosityLevel verb);

	#define	LogAll(logType, ...)	((*(Su2Log::myLog))(logType,  __FILE__, __LINE__, __VA_ARGS__))
	#define	LogDebug(...)		((*(Su2Log::myLog))(DebugLog, __FILE__, __LINE__, __VA_ARGS__))
	#define	LogError(...)		((*(Su2Log::myLog))(ErrorLog, __FILE__, __LINE__, __VA_ARGS__))
	#define	LogMsg(verb, ...)	do { if (Su2Log::myLog->Verbosity() >= verb) { ((*(Su2Log::myLog))(MsgLog, __FILE__, __LINE__, __VA_ARGS__)); } } while(0)
	#define LogOut(...) 		do { if (!Su2Comms::rank) { printf(__VA_ARGS__); fflush(stdout); } } while(0)
	#define	LogFlush()		(Su2Log::myLog->flushLog())
	#define	Timer()			(Su2Log::myLog->runTime())
#endif
