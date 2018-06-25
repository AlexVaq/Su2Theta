#include <cuda.h>
#include <cuda_runtime.h>

inline void CudaCheckError (const char *file, const int line)
{
	cudaError err = cudaDeviceSynchronize();

	if (err != cudaSuccess)
		LogError ("CudaCheckError() with sync failed: %s\n", cudaGetErrorString (err));
}
