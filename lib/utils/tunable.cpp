#ifndef	__SU2_TUNE
	#define	__SU2_TUNE

	#include <enumFields.h>
	#include <utils/utils.h>

	namespace Su2Tune {

		void    Tune (Tunable &func) {

			Su2Prof::Profiler &prof = Su2Prof::getProfiler(ProfTuner);
			prof.start();

			int  nThreads = 1;
			bool newFile  = false, found = false;

//Profiler &prof = getProfiler(PROF_TUNER);

			std::chrono::high_resolution_clock::time_point start, end;
			size_t bestTime, lastTime, cTime;

			LogMsg (VerbHigh, "Started tuner");
//prof.start();

			//if (field->Device() == DEV_CPU)
				func.ResetBlockSize();
			//else
			//	prop->InitBlockSize(field->Length(), field->Depth(), field->DataSize(), field->DataAlign(), true);

			/*      Check for a cache file  */
			//if (myRank == 0)
			{
        			FILE *cacheFile;
				char tuneName[2048];
				sprintf (tuneName, "Su2Cache.dat");
				if ((cacheFile = fopen(tuneName, "r")) == nullptr) {
					LogMsg (VerbNormal, "Missing tuning cache file %s, will create a new one", tuneName);
					newFile = true;
				} else {
				        size_t		rMpi, rThreads, rLs, rLt;
				        unsigned int	rBx, rBy, rBz, rBt, fType, sWth;//, myField = (field->Field() == FIELD_SAXION) ? 0 : 1;
				        char		mDev[8], mFunc[256], mSize[16], mCol[16];

					std::string tDev("Cpu");//field->Device() == DEV_GPU ? "Gpu" : "Cpu");

					do {
						fscanf (cacheFile, "%s %s %lu %s %s %lu %lu %u %u %u %u %u\n", &mFunc, &mCol, &sWth, &mSize, &mDev, &rMpi, &rThreads, &fType, &rBx, &rBy, &rBz, &rBt);
						std::string fDev (mDev);
						std::string nFunc(mFunc);
						nFunc += std::string(" ") + std::string(mCol) + std::string("\t") + std::to_string(sWth) + std::string(" ") + std::string(mSize);

						if (nFunc == func.Name() && rThreads == omp_get_max_threads()) { // && rMpi == commSize() && fDev == tDev) { // && rLs == func.SLength() && rLt == func.TLength()) {
							if (rBx <= func.MaxBlockX() && rBy <= func.MaxBlockY() && rBz <= func.MaxBlockZ() && rBt <= func.MaxBlockT()) {
								found = true;
								func.SetBlockX(rBx);
								func.SetBlockY(rBy);
								func.SetBlockZ(rBz);
								func.SetBlockT(rBt);
								func.UpdateBestBlock();
							}
						}
					}	while(!feof(cacheFile) && !found);
					fclose (cacheFile);
				}
			}

			//MPI_Bcast (&found, sizeof(found), MPI_BYTE, 0, MPI_COMM_WORLD);
			//commSync();

			// If a cache file was found, we broadcast the best block and exit
			if (found) {
				/*
				unsigned int block[4];

				if (myRank == 0) {
					block[0] = act->TunedBlockX();
					block[1] = act->TunedBlockY();
					block[2] = act->TunedBlockZ();
					block[3] = act->TunedBlockT();
				}

				MPI_Bcast (&block, sizeof(int)*3, MPI_BYTE, 0, MPI_COMM_WORLD);
				commSync();

				if (myRank != 0) {
					act->SetBlockX(block[0]);
					act->SetBlockY(block[1]);
					act->SetBlockZ(block[2]);
					act->SetBlockT(block[3]);
					act->UpdateBestBlock();
				}
				*/
				LogMsg (VerbNormal, "Tuned values read from cache file. Best block %u x %u x %u x %u", func.TunedBlockX(), func.TunedBlockY(), func.TunedBlockZ(), func.TunedBlockT());
				LogMsg (VerbHigh,   "Chosen block %u x %u x %u x %u", func.BlockX(), func.BlockY(), func.BlockZ(), func.BlockT());
				func.Tune();
				//prof.stop();
				//prof.add(prop->Name(), 0., 0.);

				return;
			}
			// Otherwise we start tuning

			func.SaveState();

			start = std::chrono::high_resolution_clock::now();
			func.Run();
			end   = std::chrono::high_resolution_clock::now();

			cTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
			bestTime = cTime;
			/*
			// If there is an error in GPU propagation, we set the time to an absurd value
			#ifdef USE_GPU
			if (field->Device() == DEV_GPU) {
				auto gErr = cudaGetLastError();

				if (gErr != cudaSuccess)
					cTime = std::numeric_limits<std::size_t>::max();
			}
			#endif
			*/

			//MPI_Allreduce(&cTime, &bestTime, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);

			//if (field->Device() == DEV_GPU && cTime == std::numeric_limits<std::size_t>::max())
				//LogMsg (VERB_HIGH, "Block %u x %u x %u gave an error and couldn't run on the GPU", prop->BlockX(), prop->BlockY(), prop->BlockZ());
			//else
				LogMsg (VerbHigh, "Block %u x %u x %u x %u done in %lu ns", func.BlockX(), func.BlockY(), func.BlockZ(), func.BlockT(), bestTime);
				func.AdvanceBlockSize();

			while (!func.IsTuned()) {

				start = std::chrono::high_resolution_clock::now();
				func.Run();
				end   = std::chrono::high_resolution_clock::now();

				cTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
				lastTime = cTime;
				//start = std::chrono::high_resolution_clock::now();
				//func.Run();
				//end   = std::chrono::high_resolution_clock::now();

				//cTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
				//printf("Last %zu vs Current %zu\n", lastTime, cTime);
				//MPI_Allreduce(&cTime, &lastTime, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);

				//if (field->Device() == DEV_GPU && cTime == std::numeric_limits<std::size_t>::max())
					//LogMsg (VERB_HIGH, "Block %u x %u x %u gave an error and couldn't run on the GPU", act->BlockX(), act->BlockY(), act->BlockZ());
                		//else
					LogMsg (VerbHigh, "Block %u x %u x %u x %u done in %zu ns", func.BlockX(), func.BlockY(), func.BlockZ(), func.BlockT(), lastTime);

		                if (lastTime < bestTime) {
					bestTime = lastTime;
					func.UpdateBestBlock();
					LogMsg (VerbHigh, "Best block updated");
				}

				func.AdvanceBlockSize();
			}

			func.SetBestBlock();
			LogMsg (VerbNormal, "%s tuned! Best block %u x %u x %u x %u in %zu ns", func.Name().c_str(), func.TunedBlockX(), func.TunedBlockY(), func.TunedBlockZ(), func.TunedBlockT(), bestTime);
			func.Reset();
			func.RestoreState();

			/*      Write cache file if necessary, block of rank 0 prevails         */

			//if (myRank == 0)
			{
				FILE *cacheFile;
				char tuneName[2048];
				sprintf (tuneName, "Su2Cache.dat");

				// We distinguish between opening and appending a new line
				if (!newFile) {
					if ((cacheFile = fopen(tuneName, "a")) == nullptr) {
						LogError ("Error: can't open cache file, can't save tuning results");
						//commSync();
						//prof.stop();
						//prof.add(prop->Name(), 0., 0.);
					}
				} else {
					if ((cacheFile = fopen(tuneName, "w")) == nullptr) {
						LogError ("Error: can't create cache file, can't save tuning results");
						//commSync();
						//prof.stop();
						//prof.add(prop->Name(), 0., 0.);
					}
				}

				unsigned int fType = 0; //(field->Field() == FIELD_SAXION) ? 0 : 1;
				std::string myDev("Cpu");//field->Device() == DEV_GPU ? "Gpu" : "Cpu");
				fprintf (cacheFile, "%s %s %lu %lu %u %u %u %u %u\n", func.Name().c_str(), myDev.c_str(), 1, omp_get_max_threads(),
						fType, func.TunedBlockX(), func.TunedBlockY(), func.TunedBlockZ(), func.TunedBlockT());
				fclose  (cacheFile);
			}

			prof.stop();
			std::string myName = std::string("Tuner ") + func.Name();
			prof.add(myName, 0.0, 0.0);

			LogMsg  (VerbHigh, "%s took %le ns", myName.c_str(), prof.Prof()[myName].DTime());
		}
	}
#endif

