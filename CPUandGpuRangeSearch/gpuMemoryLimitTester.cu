#include "gpuMemoryLimitTester.cuh"


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
	if (code != cudaSuccess)
	{
		std::fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}


void populateArrayWithRandomNumbers(int *array, int size)
{
	for (int i = 0; i < size; i++)
	{
		array[i] = rand();
	}
};

void checkGpuMemoryLimit()
{
	int memoryAmount =30000;
	//int memoryIncrement=1000;

	int* host_array;
	int* device_array;

	host_array = (int *)malloc(3000 * sizeof(int));
	populateArrayWithRandomNumbers(host_array, 3000);


	for (int i = 0; i < 300000; i++)
	{
		std::cout << "Attempting to allocate gpu memory size of: " << memoryAmount << std::endl;

		gpuErrchk(cudaMalloc((void**)&device_array, memoryAmount*sizeof(int)));
		//host_array = (int *)malloc(3000*sizeof(int));
		//populateArrayWithRandomNumbers(host_array, 3000);
		cudaDeviceSynchronize();
		gpuErrchk(cudaMemcpy(device_array, host_array, 3000 * sizeof(int), cudaMemcpyHostToDevice));
		cudaDeviceSynchronize();


		gpuErrchk(cudaFree(device_array));
		std::cout << "Succesfully copied data to a: " << memoryAmount << " gpu array" << std::endl << std::endl;
		memoryAmount = memoryAmount * 2;//+ memoryIncrement;
	}

};
