#include "constructKdTreesOnLoadedDataWithGPU.cuh"
__device__ kdTask taskList[20000];

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
	if (code != cudaSuccess)
	{
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}

__global__ void inversebackReferenceListToReferenceList(int *device_referenceList, int *device_backReferenceList)
{
	int i = threadIdx.x + blockDim.x * blockIdx.x;
	int position = device_backReferenceList[i];
	device_referenceList[position] = i;
}



void populate_X_KdNodeLocal(int taskNo, kdTask * pendingList, int * d_xBackPositionReferenceList, int totalNumOfConsecutiveTasks, /*int searchRange,*/ int* kdArray, int* d_xPositionReferenceList, int*d_yPositionReferenceList, int* d_zPositionReferenceList)
{
	kdTask task = pendingList[taskNo];

	if (task.kdDestination < 0)
		return;

	int activeElementCount = 0;
	int secondCount = -1;
	int elementToBeExtracted = -1;
	int posOfExtractTarget = -1;
	int currentxBackPosition;


	for (int i = task.minX; i < task.maxX; i++)
	{
		currentxBackPosition = d_xBackPositionReferenceList[i];
		if (d_xPositionReferenceList[currentxBackPosition] != -1)
		{
			if (d_yPositionReferenceList[currentxBackPosition] <= task.maxY && d_yPositionReferenceList[currentxBackPosition] >= task.minY)
			{
				if (d_zPositionReferenceList[currentxBackPosition] <= task.maxZ && d_zPositionReferenceList[currentxBackPosition] >= task.minZ)
				{
					activeElementCount++;
				}
			}
		}
	}

	//If there are no elements in the set that ahvent already been assigned a slot.
	if (activeElementCount == 0)
	{
		if (task.kdDestination>-1)
			kdArray[task.kdDestination] = -1;
		return;


	}


	if (activeElementCount % 2 == 0){
		elementToBeExtracted = activeElementCount / 2 - 1;
	}
	else{
		elementToBeExtracted = (activeElementCount - 1) / 2;
	}


	int i = task.minX;
	while (posOfExtractTarget == -1)
	{
		currentxBackPosition = d_xBackPositionReferenceList[i];
		if (d_xPositionReferenceList[currentxBackPosition] != -1)
			if (d_yPositionReferenceList[currentxBackPosition] <= task.maxY && d_yPositionReferenceList[currentxBackPosition] >= task.minY)
				if (d_zPositionReferenceList[currentxBackPosition] <= task.maxZ && d_zPositionReferenceList[currentxBackPosition] >= task.minZ)
					secondCount++;
		if (secondCount == elementToBeExtracted)
		{
			posOfExtractTarget = i;
			secondCount++;
		}
		i++;
	}




	kdTask lowerValTask;

	lowerValTask.maxX = posOfExtractTarget;
	lowerValTask.maxY = task.maxY;
	lowerValTask.maxZ = task.maxZ;
	lowerValTask.minX = task.minX;
	lowerValTask.minY = task.minY;
	lowerValTask.minZ = task.minZ;
	lowerValTask.kdDestination = task.kdDestination * 2 + 1;

	kdTask higherValTask;

	higherValTask.maxX = task.maxX;
	higherValTask.maxY = task.maxY;
	higherValTask.maxZ = task.maxZ;
	higherValTask.minY = task.minY;
	higherValTask.minX = posOfExtractTarget + 1;
	higherValTask.minY = task.minY;
	higherValTask.minZ = task.minZ;
	higherValTask.kdDestination = task.kdDestination * 2 + 2;

	pendingList[taskNo] = lowerValTask;
	pendingList[taskNo + totalNumOfConsecutiveTasks] = higherValTask;

	kdArray[task.kdDestination] = d_xBackPositionReferenceList[posOfExtractTarget];
	d_xPositionReferenceList[d_xBackPositionReferenceList[posOfExtractTarget]] = -1;
	d_yPositionReferenceList[d_xBackPositionReferenceList[posOfExtractTarget]] = -1;
	d_zPositionReferenceList[d_xBackPositionReferenceList[posOfExtractTarget]] = -1;

	//d_xBackPositionReferenceList[posOfExtractTarget] = -9999999;
}

void populate_Y_KdNodeLocal(int taskNo, kdTask * pendingList, int * d_yBackPositionReferenceList, int totalNumOfConsecutiveTasks, /*int searchRange,*/ int* kdArray, int* d_xPositionReferenceList, int*d_yPositionReferenceList, int* d_zPositionReferenceList)
{

	kdTask task = pendingList[taskNo];
	int activeElementCount = 0;
	int secondCount = -1;
	int elementToBeExtracted = -1;
	int posOfExtractTarget = -1;
	int currentyBackPosition;

	for (int i = task.minY; i < task.maxY; i++)
	{
		currentyBackPosition = d_yBackPositionReferenceList[i];
		if (d_yPositionReferenceList[currentyBackPosition] != -1)
			if (d_xPositionReferenceList[currentyBackPosition] <= task.maxX && d_xPositionReferenceList[currentyBackPosition] >= task.minX)
				if (d_zPositionReferenceList[currentyBackPosition] <= task.maxZ && d_zPositionReferenceList[currentyBackPosition] >= task.minZ)
					activeElementCount++;
	}


	if (activeElementCount == 0)
	{
		if (task.kdDestination>-1)
			kdArray[task.kdDestination] = -1;
		return;

	}


	if (activeElementCount % 2 == 0)
		elementToBeExtracted = activeElementCount / 2 - 1;
	else
		elementToBeExtracted = (activeElementCount - 1) / 2;


	int i = task.minY;
	while (posOfExtractTarget == -1)
	{
		currentyBackPosition = d_yBackPositionReferenceList[i];
		if (d_yPositionReferenceList[currentyBackPosition] != -1)
			if (d_xPositionReferenceList[currentyBackPosition] <= task.maxX && d_xPositionReferenceList[currentyBackPosition] >= task.minX)
				if (d_zPositionReferenceList[currentyBackPosition] <= task.maxZ && d_zPositionReferenceList[currentyBackPosition] >= task.minZ)
					secondCount++;
		if (secondCount == elementToBeExtracted)
		{
			posOfExtractTarget = i;
			secondCount++;
		}
		i++;
	}



	kdTask lowerValTask;

	lowerValTask.maxX = task.maxX;
	lowerValTask.maxY = posOfExtractTarget;
	lowerValTask.maxZ = task.maxZ;
	lowerValTask.minX = task.minX;
	lowerValTask.minY = task.minY;
	lowerValTask.minZ = task.minZ;
	lowerValTask.kdDestination = task.kdDestination * 2 + 1;

	kdTask higherValTask;

	higherValTask.maxX = task.maxX;
	higherValTask.maxY = task.maxY;
	higherValTask.maxZ = task.maxZ;
	higherValTask.minX = task.minX;
	higherValTask.minY = posOfExtractTarget + 1;
	higherValTask.minZ = task.minZ;
	higherValTask.kdDestination = task.kdDestination * 2 + 2;

	pendingList[taskNo] = lowerValTask;
	pendingList[taskNo + totalNumOfConsecutiveTasks] = higherValTask;

	kdArray[task.kdDestination] = d_yBackPositionReferenceList[posOfExtractTarget];
	d_xPositionReferenceList[d_yBackPositionReferenceList[posOfExtractTarget]] = -1;
	d_yPositionReferenceList[d_yBackPositionReferenceList[posOfExtractTarget]] = -1;
	d_zPositionReferenceList[d_yBackPositionReferenceList[posOfExtractTarget]] = -1;


}

void populate_Z_KdNodeLocal(int taskNo, kdTask * pendingList, int * d_zBackPositionReferenceList, int totalNumOfConsecutiveTasks, /*int searchRange,*/ int* kdArray, int* d_xPositionReferenceList, int*d_yPositionReferenceList, int* d_zPositionReferenceList)
{

	kdTask task = pendingList[taskNo];

	if (task.kdDestination < 0)
		return;

	int activeElementCount = 0;
	int secondCount = -1;
	int elementToBeExtracted = -1;
	int posOfExtractTarget = -1;
	int currentzBackPosition;

	for (int i = task.minZ; i < task.maxZ; i++)
	{
		currentzBackPosition = d_zBackPositionReferenceList[i];
		if (d_zPositionReferenceList[currentzBackPosition] != -1)
		{
			if (d_yPositionReferenceList[currentzBackPosition] <= task.maxY && d_yPositionReferenceList[currentzBackPosition] >= task.minY)
			{
				if (d_xPositionReferenceList[currentzBackPosition] <= task.maxX && d_xPositionReferenceList[currentzBackPosition] >= task.minX)
				{
					activeElementCount++;
				}
			}
		}
	}




	if (activeElementCount == 0)
	{
		if (task.kdDestination>-1)
			kdArray[task.kdDestination] = -1;
		return;


	}


	if (activeElementCount != 0)
	{
		if (activeElementCount % 2 == 0)
			elementToBeExtracted = activeElementCount / 2 - 1;
		else
			elementToBeExtracted = (activeElementCount - 1) / 2;

		int i = task.minZ;
		while (posOfExtractTarget == -1)
		{
			currentzBackPosition = d_zBackPositionReferenceList[i];
			if (d_zPositionReferenceList[currentzBackPosition] != -1)
				if (d_yPositionReferenceList[currentzBackPosition] <= task.maxY && d_yPositionReferenceList[currentzBackPosition] >= task.minY)
					if (d_xPositionReferenceList[currentzBackPosition] <= task.maxX && d_xPositionReferenceList[currentzBackPosition] >= task.minX)
						secondCount++;
			if (secondCount == elementToBeExtracted)
			{
				posOfExtractTarget = i;
				secondCount++;
			}
			i++;
		}



		kdTask lowerValTask;

		lowerValTask.maxX = task.maxX;
		lowerValTask.maxY = task.maxY;
		lowerValTask.maxZ = posOfExtractTarget;
		lowerValTask.minX = task.minX;
		lowerValTask.minY = task.minY;
		lowerValTask.minZ = task.minZ;
		lowerValTask.kdDestination = task.kdDestination * 2 + 1;

		kdTask higherValTask;

		higherValTask.maxX = task.maxX;
		higherValTask.maxY = task.maxY;
		higherValTask.maxZ = task.maxZ;
		higherValTask.minX = task.minX;
		higherValTask.minY = task.minY;
		higherValTask.minZ = posOfExtractTarget + 1;
		higherValTask.kdDestination = task.kdDestination * 2 + 2;

		pendingList[taskNo] = lowerValTask;
		pendingList[taskNo + totalNumOfConsecutiveTasks] = higherValTask;

		kdArray[task.kdDestination] = d_zBackPositionReferenceList[posOfExtractTarget];
		d_xPositionReferenceList[d_zBackPositionReferenceList[posOfExtractTarget]] = -1;
		d_yPositionReferenceList[d_zBackPositionReferenceList[posOfExtractTarget]] = -1;
		d_zPositionReferenceList[d_zBackPositionReferenceList[posOfExtractTarget]] = -1;
	}
	else
		kdArray[task.kdDestination] = -2;

}


int sizeOfFinalKdLevelLocal(int chainSize, int &numberOfLevels)
{
	numberOfLevels = 0;
	float levelSize = 0.5;
	int remainingChain = chainSize;
	while (remainingChain > 0)
	{
		numberOfLevels++;
		levelSize = levelSize * 2;
		remainingChain = remainingChain - levelSize;
	}

	return int(levelSize);
}



__global__ void bitonic_sort_step(int *device_coordArray, int j, int k, int *device_backReferenceList)
{
	unsigned int i, ixj; /* Sorting partners: i and ixj */
	i = threadIdx.x + blockDim.x * blockIdx.x;
	ixj = i^j;

	i = i;
	ixj = ixj;

	/* The threads with the lowest ids sort the array. */
	if ((ixj) > i) {
		if ((i&k) == 0) {
			/* Sort ascending */
			if (device_coordArray[i] > device_coordArray[ixj]) {
				/* exchange(i,ixj); */
				int temp = device_coordArray[i];
				device_coordArray[i] = device_coordArray[ixj];
				device_coordArray[ixj] = temp;
				int temptwo = device_backReferenceList[i];
				device_backReferenceList[i] = device_backReferenceList[ixj];
				device_backReferenceList[ixj] = temptwo;
			}
		}
		if ((i&k) != 0) {
			/* Sort descending */
			if (device_coordArray[i] < device_coordArray[ixj]) {
				/* exchange(i,ixj); */
				int temp = device_coordArray[i];
				device_coordArray[i] = device_coordArray[ixj];
				device_coordArray[ixj] = temp;
				temp = device_backReferenceList[i];
				device_backReferenceList[i] = device_backReferenceList[ixj];
				device_backReferenceList[ixj] = temp;
			}
		}
	}
}


void constructKdTreesOnLoadedDataOnCPUWithGPUSorting(ProteinDataHandler &ProteinData)
{
	dataLoaderWithGpuExtensions kdTreeConstructionArtifacts;
	cudaStream_t streams[3];
	cudaStreamCreate(&streams[0]);
	cudaStreamCreate(&streams[1]);
	cudaStreamCreate(&streams[2]);
	kdTreeConstructionArtifacts.xValues = new int[16384 * 10];
	kdTreeConstructionArtifacts.yValues = new int[16384 * 10];
	kdTreeConstructionArtifacts.zValues = new int[16384 * 10];
	kdTreeConstructionArtifacts.xBackPositionReferenceList = new int[16384 * 10];
	kdTreeConstructionArtifacts.yBackPositionReferenceList = new int[16384 * 10];
	kdTreeConstructionArtifacts.zBackPositionReferenceList = new int[16384 * 10];
	kdTreeConstructionArtifacts.xForwardPositionReferenceList = new int[16384 * 10];
	kdTreeConstructionArtifacts.yForwardPositionReferenceList = new int[16384 * 10];
	kdTreeConstructionArtifacts.zForwardPositionReferenceList = new int[16384 * 10];
	kdTreeConstructionArtifacts.kdTreeholder = new int[16384 * 10 * 2];
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_xBackPositionReferenceList, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_yBackPositionReferenceList, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_zBackPositionReferenceList, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_xForwardPositionReferenceList, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_yForwardPositionReferenceList, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_zForwardPositionReferenceList, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_xCoordsHolder, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_yCoordsHolder, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_zCoordsHolder, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_xForwardPositionReferenceList, sizeof(int) * 16384 * 10));
	//gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_kdTree, sizeof(int) * 16384 * 10));

	for (int i = 0; i<5; i++)
	{
		for (int j = 0; j<ProteinData.ProteinDataHolder[i].heldEntries; j++)
		{
			//Load data of next protein and reset intermediate value arrays
			for (int y = 0; y<ProteinData.ProteinDataHolder[i].MaxEntrySize; y++)
			{
				kdTreeConstructionArtifacts.xValues[y] = ProteinData.ProteinDataHolder[i].xCoordsSets[y + ProteinData.ProteinDataHolder[i].MaxEntrySize*j];
				kdTreeConstructionArtifacts.yValues[y] = ProteinData.ProteinDataHolder[i].yCoordsSets[y + ProteinData.ProteinDataHolder[i].MaxEntrySize*j];
				kdTreeConstructionArtifacts.zValues[y] = ProteinData.ProteinDataHolder[i].zCoordsSets[y + ProteinData.ProteinDataHolder[i].MaxEntrySize*j];
				kdTreeConstructionArtifacts.xBackPositionReferenceList[y] = y;
				kdTreeConstructionArtifacts.yBackPositionReferenceList[y] = y;
				kdTreeConstructionArtifacts.zBackPositionReferenceList[y] = y;
				kdTreeConstructionArtifacts.xForwardPositionReferenceList[y] = y;
				kdTreeConstructionArtifacts.yForwardPositionReferenceList[y] = y;
				kdTreeConstructionArtifacts.zForwardPositionReferenceList[y] = y;
			}
			std::fill(kdTreeConstructionArtifacts.kdTreeholder, kdTreeConstructionArtifacts.kdTreeholder + ProteinData.ProteinDataHolder[i].MaxEntrySize * 2, -1);
			resetReverseAndForwardArrays << <1, 1, 0, streams[1] >> >(kdTreeConstructionArtifacts.d_xBackPositionReferenceList, kdTreeConstructionArtifacts.d_xForwardPositionReferenceList, kdTreeConstructionArtifacts.d_yBackPositionReferenceList, kdTreeConstructionArtifacts.d_yForwardPositionReferenceList, kdTreeConstructionArtifacts.d_zBackPositionReferenceList, kdTreeConstructionArtifacts.d_zForwardPositionReferenceList, 16384 * 10);

			int threads = 512;
			int blocks = ProteinData.ProteinDataHolder[i].MaxEntrySize / 512;


			gpuErrchk(cudaMemsetAsync(kdTreeConstructionArtifacts.d_xCoordsHolder, 999999, (ProteinData.ProteinDataHolder[i].proteinLengthCounts[j] - 1), streams[0]));
			gpuErrchk(cudaMemsetAsync(kdTreeConstructionArtifacts.d_yCoordsHolder, 999999, (ProteinData.ProteinDataHolder[i].proteinLengthCounts[j] - 1), streams[1]));
			gpuErrchk(cudaMemsetAsync(kdTreeConstructionArtifacts.d_zCoordsHolder, 999999, (ProteinData.ProteinDataHolder[i].proteinLengthCounts[j] - 1), streams[2]));
			Parallel_Device_Side_bitonic_sort(kdTreeConstructionArtifacts.xValues, kdTreeConstructionArtifacts.xBackPositionReferenceList, kdTreeConstructionArtifacts.d_xForwardPositionReferenceList, kdTreeConstructionArtifacts.d_xBackPositionReferenceList, kdTreeConstructionArtifacts.d_xCoordsHolder, kdTreeConstructionArtifacts.yValues, kdTreeConstructionArtifacts.yBackPositionReferenceList, kdTreeConstructionArtifacts.d_yForwardPositionReferenceList, kdTreeConstructionArtifacts.d_yBackPositionReferenceList, kdTreeConstructionArtifacts.d_yCoordsHolder, kdTreeConstructionArtifacts.zValues, kdTreeConstructionArtifacts.zBackPositionReferenceList, kdTreeConstructionArtifacts.d_zForwardPositionReferenceList, kdTreeConstructionArtifacts.d_zBackPositionReferenceList, kdTreeConstructionArtifacts.d_zCoordsHolder, blocks, threads);



			//sortReferenceArrayByQuickSort(kdTreeConstructionArtifacts.xValues, kdTreeConstructionArtifacts.xBackPositionReferenceList, 0, (ProteinData.ProteinDataHolder[i].proteinLengthCounts[j] - 1));
			//setforwardReferenceArrayFromBack(kdTreeConstructionArtifacts.xForwardPositionReferenceList, kdTreeConstructionArtifacts.xBackPositionReferenceList, ProteinData.ProteinDataHolder[i].proteinLengthCounts[j]);

			//sortReferenceArrayByQuickSort(kdTreeConstructionArtifacts.yValues, kdTreeConstructionArtifacts.yBackPositionReferenceList, 0, (ProteinData.ProteinDataHolder[i].proteinLengthCounts[j] - 1));
			//setforwardReferenceArrayFromBack(kdTreeConstructionArtifacts.yForwardPositionReferenceList, kdTreeConstructionArtifacts.yBackPositionReferenceList, ProteinData.ProteinDataHolder[i].proteinLengthCounts[j]);

			//sortReferenceArrayByQuickSort(kdTreeConstructionArtifacts.zValues, kdTreeConstructionArtifacts.zBackPositionReferenceList, 0, (ProteinData.ProteinDataHolder[i].proteinLengthCounts[j] - 1));
			//setforwardReferenceArrayFromBack(kdTreeConstructionArtifacts.zForwardPositionReferenceList, kdTreeConstructionArtifacts.zBackPositionReferenceList, ProteinData.ProteinDataHolder[i].proteinLengthCounts[j]);

			//Need to get values back from the gpu

			//construct kdTree
			size_t coordinatesSize = blocks* threads*sizeof(int);

			gpuErrchk(cudaMemcpyAsync(kdTreeConstructionArtifacts.xForwardPositionReferenceList, kdTreeConstructionArtifacts.d_xForwardPositionReferenceList, coordinatesSize, cudaMemcpyDeviceToHost, streams[0]));
			gpuErrchk(cudaMemcpyAsync(kdTreeConstructionArtifacts.yForwardPositionReferenceList, kdTreeConstructionArtifacts.d_yForwardPositionReferenceList, coordinatesSize, cudaMemcpyDeviceToHost, streams[1]));
			gpuErrchk(cudaMemcpyAsync(kdTreeConstructionArtifacts.zForwardPositionReferenceList, kdTreeConstructionArtifacts.d_zForwardPositionReferenceList, coordinatesSize, cudaMemcpyDeviceToHost, streams[2]));

			cudaDeviceSynchronize();

			kdTree_constructorLocal(kdTreeConstructionArtifacts.kdTreeholder, kdTreeConstructionArtifacts.xForwardPositionReferenceList, kdTreeConstructionArtifacts.xBackPositionReferenceList, kdTreeConstructionArtifacts.yForwardPositionReferenceList, kdTreeConstructionArtifacts.yBackPositionReferenceList, kdTreeConstructionArtifacts.zForwardPositionReferenceList, kdTreeConstructionArtifacts.zBackPositionReferenceList, ProteinData.ProteinDataHolder[i].proteinLengthCounts[j], ProteinData.ProteinDataHolder[i].MaxEntrySize);
	
			
			
			//int X = 4;

			//move kdTree to main storage arrays
			for (int y = 0; y<ProteinData.ProteinDataHolder[i].MaxEntrySize * 2; y++)
			{
				ProteinData.ProteinDataHolder[i].kdTrees[y + j* ProteinData.ProteinDataHolder[i].KdTreeSize] = kdTreeConstructionArtifacts.kdTreeholder[y];

			}
		}
	}





}


void constructKdTreesOnLoadedDataOnGPU(ProteinDataHandler &ProteinData)
{
	dataLoaderWithGpuExtensions kdTreeConstructionArtifacts;
	cudaStream_t streams[3];
	cudaStreamCreate(&streams[0]);
	cudaStreamCreate(&streams[1]);
	cudaStreamCreate(&streams[2]);
	kdTreeConstructionArtifacts.xValues = new int[16384 * 10];
	kdTreeConstructionArtifacts.yValues = new int[16384 * 10];
	kdTreeConstructionArtifacts.zValues = new int[16384 * 10];
	kdTreeConstructionArtifacts.xBackPositionReferenceList = new int[16384 * 10];
	kdTreeConstructionArtifacts.yBackPositionReferenceList = new int[16384 * 10];
	kdTreeConstructionArtifacts.zBackPositionReferenceList = new int[16384 * 10];
	kdTreeConstructionArtifacts.xForwardPositionReferenceList = new int[16384 * 10];
	kdTreeConstructionArtifacts.yForwardPositionReferenceList = new int[16384 * 10];
	kdTreeConstructionArtifacts.zForwardPositionReferenceList = new int[16384 * 10];
	kdTreeConstructionArtifacts.kdTreeholder = new int[16384 * 10 * 2];
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_xBackPositionReferenceList, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_yBackPositionReferenceList, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_zBackPositionReferenceList, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_xForwardPositionReferenceList, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_yForwardPositionReferenceList, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_zForwardPositionReferenceList, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_xCoordsHolder, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_yCoordsHolder, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_zCoordsHolder, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_xForwardPositionReferenceList, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_kdTree, sizeof(int) * 16384 * 10));

	for (int i = 0; i<5; i++)
	{
		for (int j = 0; j<ProteinData.ProteinDataHolder[i].heldEntries; j++)
		{
std::cout<<"Processing entry number: "<<j<<" "<<std::endl;
			//Load data of next protein and reset intermediate value arrays
			for (int y = 0; y<ProteinData.ProteinDataHolder[i].MaxEntrySize; y++)
			{
				kdTreeConstructionArtifacts.xValues[y] = ProteinData.ProteinDataHolder[i].xCoordsSets[y + ProteinData.ProteinDataHolder[i].MaxEntrySize*j];
				kdTreeConstructionArtifacts.yValues[y] = ProteinData.ProteinDataHolder[i].yCoordsSets[y + ProteinData.ProteinDataHolder[i].MaxEntrySize*j];
				kdTreeConstructionArtifacts.zValues[y] = ProteinData.ProteinDataHolder[i].zCoordsSets[y + ProteinData.ProteinDataHolder[i].MaxEntrySize*j];
				kdTreeConstructionArtifacts.xBackPositionReferenceList[y] = y;
				kdTreeConstructionArtifacts.yBackPositionReferenceList[y] = y;
				kdTreeConstructionArtifacts.zBackPositionReferenceList[y] = y;
				kdTreeConstructionArtifacts.xForwardPositionReferenceList[y] = y;
				kdTreeConstructionArtifacts.yForwardPositionReferenceList[y] = y;
				kdTreeConstructionArtifacts.zForwardPositionReferenceList[y] = y;
			}
			std::fill(kdTreeConstructionArtifacts.kdTreeholder, kdTreeConstructionArtifacts.kdTreeholder + ProteinData.ProteinDataHolder[i].MaxEntrySize * 2, -1);
			resetReverseAndForwardArrays << <1, 1, 0, streams[1] >> >(kdTreeConstructionArtifacts.d_xBackPositionReferenceList, kdTreeConstructionArtifacts.d_xForwardPositionReferenceList, kdTreeConstructionArtifacts.d_yBackPositionReferenceList, kdTreeConstructionArtifacts.d_yForwardPositionReferenceList, kdTreeConstructionArtifacts.d_zBackPositionReferenceList, kdTreeConstructionArtifacts.d_zForwardPositionReferenceList, 16384 * 10);

			int threads = 512;
			int blocks = ProteinData.ProteinDataHolder[i].MaxEntrySize / 512;

			//std::cout << "Working on pair:" << i << " " << j << std::endl;
			gpuErrchk(cudaMemsetAsync(kdTreeConstructionArtifacts.d_xCoordsHolder, 999999, (ProteinData.ProteinDataHolder[i].proteinLengthCounts[j] - 1), streams[0]));
			gpuErrchk(cudaMemsetAsync(kdTreeConstructionArtifacts.d_yCoordsHolder, 999999, (ProteinData.ProteinDataHolder[i].proteinLengthCounts[j] - 1), streams[1]));
			gpuErrchk(cudaMemsetAsync(kdTreeConstructionArtifacts.d_zCoordsHolder, 999999, (ProteinData.ProteinDataHolder[i].proteinLengthCounts[j] - 1), streams[2]));
			Parallel_Device_Side_bitonic_sort(kdTreeConstructionArtifacts.xValues, kdTreeConstructionArtifacts.xBackPositionReferenceList, kdTreeConstructionArtifacts.d_xForwardPositionReferenceList, kdTreeConstructionArtifacts.d_xBackPositionReferenceList, kdTreeConstructionArtifacts.d_xCoordsHolder, kdTreeConstructionArtifacts.yValues, kdTreeConstructionArtifacts.yBackPositionReferenceList, kdTreeConstructionArtifacts.d_yForwardPositionReferenceList, kdTreeConstructionArtifacts.d_yBackPositionReferenceList, kdTreeConstructionArtifacts.d_yCoordsHolder, kdTreeConstructionArtifacts.zValues, kdTreeConstructionArtifacts.zBackPositionReferenceList, kdTreeConstructionArtifacts.d_zForwardPositionReferenceList, kdTreeConstructionArtifacts.d_zBackPositionReferenceList, kdTreeConstructionArtifacts.d_zCoordsHolder, blocks, threads);


			int* numOfAtoms = (int*)malloc(sizeof(int));
			numOfAtoms[0]= ProteinData.ProteinDataHolder[i].proteinLengthCounts[j];
			Device_Side_kdTree_constructor(kdTreeConstructionArtifacts.d_kdTree, kdTreeConstructionArtifacts.d_xForwardPositionReferenceList, kdTreeConstructionArtifacts.d_xBackPositionReferenceList, kdTreeConstructionArtifacts.d_yForwardPositionReferenceList, kdTreeConstructionArtifacts.d_yBackPositionReferenceList, kdTreeConstructionArtifacts.d_zForwardPositionReferenceList, kdTreeConstructionArtifacts.d_yBackPositionReferenceList, numOfAtoms, ProteinData.ProteinDataHolder[i].MaxEntrySize);
			//Device_Side_kdTree_constructor(int* d_kdArray,                       int* d_xPositionReferenceList, int* d_xBackPositionReferenceList, int* d_yPositionReferenceList, int* d_yBackPositionReferenceList, int*  d_zPositionReferenceList, int*  d_zBackPositionReferenceList, int  *numOfRealInputElements, int lengthOfHolderArrays)

			cudaDeviceSynchronize();
			gpuErrchk(cudaMemcpy(kdTreeConstructionArtifacts.kdTreeholder, kdTreeConstructionArtifacts.d_kdTree, 400 /*ProteinData.ProteinDataHolder[i].KdTreeSize*sizeof(int)*/, cudaMemcpyDeviceToHost));
		


			//kdTree_constructorLocal(kdTreeConstructionArtifacts.kdTreeholder, kdTreeConstructionArtifacts.xForwardPositionReferenceList, kdTreeConstructionArtifacts.xBackPositionReferenceList, kdTreeConstructionArtifacts.yForwardPositionReferenceList, kdTreeConstructionArtifacts.yBackPositionReferenceList, kdTreeConstructionArtifacts.zForwardPositionReferenceList, kdTreeConstructionArtifacts.zBackPositionReferenceList, ProteinData.ProteinDataHolder[i].proteinLengthCounts[j], ProteinData.ProteinDataHolder[i].MaxEntrySize);
			//int X = 4;

			//move kdTree to main storage arrays
			for (int y = 0; y<ProteinData.ProteinDataHolder[i].MaxEntrySize * 2; y++)
			{
				ProteinData.ProteinDataHolder[i].kdTrees[y + j* ProteinData.ProteinDataHolder[i].KdTreeSize] = kdTreeConstructionArtifacts.kdTreeholder[y];

			}
		}
	}







}


void Parallel_Device_Side_bitonic_sort(int *h_coordArray1, int * h_reverseValuePositionArray, int* d_forwardPositionArray, int*d_reversePositionArray, int * d_coordsHolder, int *h_coordArray2, int * h_reverseValuePositionArray2, int* d_forwardPositionArray2, int*d_reversePositionArray2, int * d_coordsHolder2, int *h_coordArray3, int * h_reverseValuePositionArray3, int* d_forwardPositionArray3, int*d_reversePositionArray3, int * d_coordsHolder3, int numOfBlock, int numOfThreads)
{


	cudaStream_t streams[3];


	cudaStreamCreate(&streams[0]);
	cudaStreamCreate(&streams[1]);
	cudaStreamCreate(&streams[2]);


	size_t coordinatesSize = (numOfBlock *numOfThreads) * sizeof(int);//Used for determining the size of the short based coord arrays for use on the gpu;


	gpuErrchk(cudaMemcpyAsync(d_coordsHolder, h_coordArray1, coordinatesSize, cudaMemcpyHostToDevice, streams[0]));
	gpuErrchk(cudaMemcpyAsync(d_coordsHolder2, h_coordArray2, coordinatesSize, cudaMemcpyHostToDevice, streams[1]));
	gpuErrchk(cudaMemcpyAsync(d_coordsHolder3, h_coordArray3, coordinatesSize, cudaMemcpyHostToDevice, streams[2]));



	gpuErrchk(cudaMemcpyAsync(d_reversePositionArray, h_reverseValuePositionArray, coordinatesSize, cudaMemcpyHostToDevice, streams[0]));
	gpuErrchk(cudaMemcpyAsync(d_reversePositionArray2, h_reverseValuePositionArray2, coordinatesSize, cudaMemcpyHostToDevice, streams[1]));
	gpuErrchk(cudaMemcpyAsync(d_reversePositionArray3, h_reverseValuePositionArray3, coordinatesSize, cudaMemcpyHostToDevice, streams[2]));


	dim3 blocks(numOfBlock, 1);
	dim3 threads(numOfThreads, 1);


	int j, k;
	// Major step 
	for (k = 2; k <= (numOfBlock *numOfThreads); k <<= 1) {
		// Minor step 
		for (j = k >> 1; j > 0; j = j >> 1) {

			bitonic_sort_step << <blocks, threads, 0, streams[0] >> >(d_coordsHolder, j, k, d_reversePositionArray); //this should leave the device_coordArray ordered (and discardable) and the device_backReferenceList pointing to where each final value location came from
			bitonic_sort_step << <blocks, threads, 0, streams[1] >> >(d_coordsHolder2, j, k, d_reversePositionArray2); //this should leave the device_coordArray ordered (and discardable) and the device_backReferenceList pointing to where each final value location came from
			bitonic_sort_step << <blocks, threads, 0, streams[2] >> >(d_coordsHolder3, j, k, d_reversePositionArray3); //this should leave the device_coordArray ordered (and discardable) and the device_backReferenceList pointing to where each final value location came from

		}
	}




	inversebackReferenceListToReferenceList << <numOfBlock, numOfThreads, 0, streams[0] >> >(d_forwardPositionArray, d_reversePositionArray);
	inversebackReferenceListToReferenceList << <numOfBlock, numOfThreads, 0, streams[1] >> >(d_forwardPositionArray2, d_reversePositionArray2);
	inversebackReferenceListToReferenceList << <numOfBlock, numOfThreads, 0, streams[2] >> >(d_forwardPositionArray3, d_reversePositionArray3);

	//cudaDeviceSynchronize();
	gpuErrchk(cudaMemcpyAsync(h_reverseValuePositionArray, d_reversePositionArray, coordinatesSize, cudaMemcpyDeviceToHost, streams[0]));
	gpuErrchk(cudaMemcpyAsync(h_reverseValuePositionArray2, d_reversePositionArray2, coordinatesSize, cudaMemcpyDeviceToHost, streams[1]));
	gpuErrchk(cudaMemcpyAsync(h_reverseValuePositionArray3, d_reversePositionArray3, coordinatesSize, cudaMemcpyDeviceToHost, streams[2]));


	gpuErrchk(cudaMemcpyAsync(h_coordArray1, d_coordsHolder, coordinatesSize, cudaMemcpyDeviceToHost, streams[0]));
	gpuErrchk(cudaMemcpyAsync(h_coordArray2, d_coordsHolder2, coordinatesSize, cudaMemcpyDeviceToHost, streams[1]));
	gpuErrchk(cudaMemcpyAsync(h_coordArray3, d_coordsHolder3, coordinatesSize, cudaMemcpyDeviceToHost, streams[2]));



	cudaStreamDestroy(streams[0]);
	cudaStreamDestroy(streams[1]);
	cudaStreamDestroy(streams[2]);

}

__global__ void resetReverseAndForwardArrays(int* d_xbackArray, int* d_xforwardArray, int* d_ybackArray, int* d_yforwardArray, int* d_zbackArray, int* d_zforwardArray, int length)
{
	for (int i = 0; i < length; i++)
	{
		d_xbackArray[i] = i;
		d_xforwardArray[i] = i;
		d_ybackArray[i] = i;
		d_yforwardArray[i] = i;
		d_zbackArray[i] = i;
		d_zforwardArray[i] = i;
	}
}


void kdTree_constructorLocal(int* kdArray, int* xPositionReferenceList, int* xBackPositionReferenceList, int* yPositionReferenceList, int* yBackPositionReferenceList, int*  zPositionReferenceList, int*  zBackPositionReferenceList, int  numOfRealInputElements, int lengthOfHolderArrays)
{
	int numOfLevels;
	sizeOfFinalKdLevelLocal(numOfRealInputElements, numOfLevels);
	int taskArraySize = lengthOfHolderArrays;
	kdTask *taskList = (kdTask*)malloc(sizeof(kdTask)*(taskArraySize * 2));

	int numberOfCoords = numOfRealInputElements;
	for (int i = 0; i < taskArraySize * 2; i++)
		taskList[i].kdDestination = -1;

	kdTask initialTask;
	initialTask.kdDestination = 0;
	initialTask.maxX = numberOfCoords;
	initialTask.maxY = numberOfCoords;
	initialTask.maxZ = numberOfCoords;
	initialTask.minX = 0;
	initialTask.minY = 0;
	initialTask.minZ = 0;

	taskList[0] = initialTask;

	int activeTasks = 1;

	//How do we confirm we processing all tasks?
	for (int i = 0; i < (numOfLevels); i++)
	{

		if (i % 3 == 0)
			for (int j = 0; j < activeTasks; j++)
				populate_X_KdNodeLocal(j, taskList, xBackPositionReferenceList, activeTasks, kdArray, xPositionReferenceList, yPositionReferenceList, zPositionReferenceList);
		else if (i % 3 == 1)
			for (int j = 0; j < activeTasks; j++)
				populate_Y_KdNodeLocal(j, taskList, yBackPositionReferenceList, activeTasks, kdArray, xPositionReferenceList, yPositionReferenceList, zPositionReferenceList);
		else if (i % 3 == 2)
			for (int j = 0; j < activeTasks; j++)
				populate_Z_KdNodeLocal(j, taskList, zBackPositionReferenceList, activeTasks, kdArray, xPositionReferenceList, yPositionReferenceList, zPositionReferenceList);
		activeTasks = activeTasks * 2;
	}
}



void Device_Side_kdTree_constructor(int* d_kdArray, int* d_xPositionReferenceList, int* d_xBackPositionReferenceList, int* d_yPositionReferenceList, int* d_yBackPositionReferenceList, int*  d_zPositionReferenceList, int*  d_zBackPositionReferenceList, int  *numOfRealInputElements, int lengthOfHolderArrays)
{
	int chainSize = lengthOfHolderArrays;
	int numberOfLevels;
	int val1 = sizeOfFinalKdLevelLocal(chainSize, numberOfLevels);

	int numberOfCoords = numOfRealInputElements[0];//number of atoms in the chain


	online_kdConstructor_master_node << <1, 1 >> >(d_kdArray, d_xPositionReferenceList, d_xBackPositionReferenceList, d_yPositionReferenceList, d_yBackPositionReferenceList, d_zPositionReferenceList, d_zBackPositionReferenceList, numberOfCoords);
	int threadcount = 1;
	int blockcount = 1;
	int activeTasks = 1;
	for (int i = 0; i < numberOfLevels - 1; i++)
	{
		if (i % 3 == 0)
			online_kdConstructor_slave_nodeX << <  blockcount, threadcount >> >(d_xBackPositionReferenceList, activeTasks, d_kdArray, d_xPositionReferenceList, d_yPositionReferenceList, d_zPositionReferenceList);
		else if (i % 3 == 1)
			online_kdConstructor_slave_nodeY << <  blockcount, threadcount >> >(d_yBackPositionReferenceList, activeTasks, d_kdArray, d_xPositionReferenceList, d_yPositionReferenceList, d_zPositionReferenceList);
		else if (i % 3 == 2)
			online_kdConstructor_slave_nodeZ << <  blockcount, threadcount >> >(d_zBackPositionReferenceList, activeTasks, d_kdArray, d_xPositionReferenceList, d_yPositionReferenceList, d_zPositionReferenceList);

		if (threadcount < 512)
			threadcount = threadcount * 2;
		else blockcount = blockcount * 2;
		activeTasks = activeTasks * 2;
	}
}





void GPUTreeWithCPULoad(ProteinDataHandler &ProteinData)
{
	dataLoaderWithGpuExtensions kdTreeConstructionArtifacts;


	kdTreeConstructionArtifacts.xValues = new int[16384 * 10];
	kdTreeConstructionArtifacts.yValues = new int[16384 * 10];
	kdTreeConstructionArtifacts.zValues = new int[16384 * 10];
	kdTreeConstructionArtifacts.xBackPositionReferenceList = new int[16384 * 10];
	kdTreeConstructionArtifacts.yBackPositionReferenceList = new int[16384 * 10];
	kdTreeConstructionArtifacts.zBackPositionReferenceList = new int[16384 * 10];
	kdTreeConstructionArtifacts.xForwardPositionReferenceList = new int[16384 * 10];
	kdTreeConstructionArtifacts.yForwardPositionReferenceList = new int[16384 * 10];
	kdTreeConstructionArtifacts.zForwardPositionReferenceList = new int[16384 * 10];
	kdTreeConstructionArtifacts.kdTreeholder = new int[16384 * 10 * 2];
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_xBackPositionReferenceList, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_yBackPositionReferenceList, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_zBackPositionReferenceList, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_xForwardPositionReferenceList, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_yForwardPositionReferenceList, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_zForwardPositionReferenceList, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_xCoordsHolder, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_yCoordsHolder, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_zCoordsHolder, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_xForwardPositionReferenceList, sizeof(int) * 16384 * 10));
	gpuErrchk(cudaMalloc((void**)&kdTreeConstructionArtifacts.d_kdTree, sizeof(int) * 16384 * 10));


	for (int i = 0; i<5; i++)
	{

		for (int j = 0; j<ProteinData.ProteinDataHolder[i].heldEntries; j++)
		{
			//Load data of next protein and reset intermediate value arrays
			for (int y = 0; y<ProteinData.ProteinDataHolder[i].MaxEntrySize; y++)
			{
				kdTreeConstructionArtifacts.xValues[y] = ProteinData.ProteinDataHolder[i].xCoordsSets[y + ProteinData.ProteinDataHolder[i].MaxEntrySize*j];
				kdTreeConstructionArtifacts.yValues[y] = ProteinData.ProteinDataHolder[i].yCoordsSets[y + ProteinData.ProteinDataHolder[i].MaxEntrySize*j];
				kdTreeConstructionArtifacts.zValues[y] = ProteinData.ProteinDataHolder[i].zCoordsSets[y + ProteinData.ProteinDataHolder[i].MaxEntrySize*j];

				kdTreeConstructionArtifacts.xBackPositionReferenceList[y] = y;
				kdTreeConstructionArtifacts.yBackPositionReferenceList[y] = y;
				kdTreeConstructionArtifacts.zBackPositionReferenceList[y] = y;
				kdTreeConstructionArtifacts.xForwardPositionReferenceList[y] = y;
				kdTreeConstructionArtifacts.yForwardPositionReferenceList[y] = y;
				kdTreeConstructionArtifacts.zForwardPositionReferenceList[y] = y;
			}
			std::fill(kdTreeConstructionArtifacts.kdTreeholder, kdTreeConstructionArtifacts.kdTreeholder + ProteinData.ProteinDataHolder[i].MaxEntrySize * 2, -1);



			sortReferenceArrayByQuickSort(kdTreeConstructionArtifacts.xValues, kdTreeConstructionArtifacts.xBackPositionReferenceList, 0, (ProteinData.ProteinDataHolder[i].proteinLengthCounts[j] - 1));
			setforwardReferenceArrayFromBack(kdTreeConstructionArtifacts.xForwardPositionReferenceList, kdTreeConstructionArtifacts.xBackPositionReferenceList, ProteinData.ProteinDataHolder[i].proteinLengthCounts[j]);

			sortReferenceArrayByQuickSort(kdTreeConstructionArtifacts.yValues, kdTreeConstructionArtifacts.yBackPositionReferenceList, 0, (ProteinData.ProteinDataHolder[i].proteinLengthCounts[j] - 1));
			setforwardReferenceArrayFromBack(kdTreeConstructionArtifacts.yForwardPositionReferenceList, kdTreeConstructionArtifacts.yBackPositionReferenceList, ProteinData.ProteinDataHolder[i].proteinLengthCounts[j]);

			sortReferenceArrayByQuickSort(kdTreeConstructionArtifacts.zValues, kdTreeConstructionArtifacts.zBackPositionReferenceList, 0, (ProteinData.ProteinDataHolder[i].proteinLengthCounts[j] - 1));
			setforwardReferenceArrayFromBack(kdTreeConstructionArtifacts.zForwardPositionReferenceList, kdTreeConstructionArtifacts.zBackPositionReferenceList, ProteinData.ProteinDataHolder[i].proteinLengthCounts[j]);




			int* numOfAtoms = (int*)malloc(sizeof(int));
			numOfAtoms[0] = ProteinData.ProteinDataHolder[i].proteinLengthCounts[j];
			
			gpuErrchk(cudaMemcpy(kdTreeConstructionArtifacts.d_xForwardPositionReferenceList, kdTreeConstructionArtifacts.xForwardPositionReferenceList, ProteinData.ProteinDataHolder[i].MaxEntrySize*sizeof(int), cudaMemcpyHostToDevice));
			gpuErrchk(cudaMemcpy(kdTreeConstructionArtifacts.d_yForwardPositionReferenceList, kdTreeConstructionArtifacts.yForwardPositionReferenceList, ProteinData.ProteinDataHolder[i].MaxEntrySize*sizeof(int), cudaMemcpyHostToDevice));
			gpuErrchk(cudaMemcpy(kdTreeConstructionArtifacts.d_zForwardPositionReferenceList, kdTreeConstructionArtifacts.zForwardPositionReferenceList, ProteinData.ProteinDataHolder[i].MaxEntrySize*sizeof(int), cudaMemcpyHostToDevice));
			gpuErrchk(cudaMemcpy(kdTreeConstructionArtifacts.d_xBackPositionReferenceList, kdTreeConstructionArtifacts.xBackPositionReferenceList, ProteinData.ProteinDataHolder[i].MaxEntrySize*sizeof(int), cudaMemcpyHostToDevice));
			gpuErrchk(cudaMemcpy(kdTreeConstructionArtifacts.d_yBackPositionReferenceList, kdTreeConstructionArtifacts.yBackPositionReferenceList, ProteinData.ProteinDataHolder[i].MaxEntrySize*sizeof(int), cudaMemcpyHostToDevice));
			gpuErrchk(cudaMemcpy(kdTreeConstructionArtifacts.d_zBackPositionReferenceList, kdTreeConstructionArtifacts.zBackPositionReferenceList, ProteinData.ProteinDataHolder[i].MaxEntrySize*sizeof(int), cudaMemcpyHostToDevice));









			Device_Side_kdTree_constructor(kdTreeConstructionArtifacts.d_kdTree, kdTreeConstructionArtifacts.d_xForwardPositionReferenceList, kdTreeConstructionArtifacts.d_xBackPositionReferenceList, kdTreeConstructionArtifacts.d_yForwardPositionReferenceList, kdTreeConstructionArtifacts.d_yBackPositionReferenceList, kdTreeConstructionArtifacts.d_zForwardPositionReferenceList, kdTreeConstructionArtifacts.d_yBackPositionReferenceList, numOfAtoms, ProteinData.ProteinDataHolder[i].MaxEntrySize);

			gpuErrchk(cudaMemcpy(kdTreeConstructionArtifacts.kdTreeholder, kdTreeConstructionArtifacts.d_kdTree, ProteinData.ProteinDataHolder[i].KdTreeSize*sizeof(int), cudaMemcpyDeviceToHost));
			//construct kdTree
			//kdTree_constructor(kdTreeConstructionArtifacts.kdTreeholder, kdTreeConstructionArtifacts.xForwardPositionReferenceList, kdTreeConstructionArtifacts.xBackPositionReferenceList, kdTreeConstructionArtifacts.yForwardPositionReferenceList, kdTreeConstructionArtifacts.yBackPositionReferenceList, kdTreeConstructionArtifacts.zForwardPositionReferenceList, kdTreeConstructionArtifacts.zBackPositionReferenceList, ProteinData.ProteinDataHolder[i].proteinLengthCounts[j], ProteinData.ProteinDataHolder[i].MaxEntrySize);


			//move kdTree to main storage arrays
			for (int y = 0; y<ProteinData.ProteinDataHolder[i].MaxEntrySize * 2; y++)
			{
				ProteinData.ProteinDataHolder[i].kdTrees[y + j* ProteinData.ProteinDataHolder[i].KdTreeSize] = kdTreeConstructionArtifacts.kdTreeholder[y];

			}

		}
	}


}





__global__ void online_kdConstructor_master_node(int * d_kdArray, int * d_xPositionReferenceList, int * d_xBackPositionReferenceList, int * d_yPositionReferenceList, int * d_yBackPositionReferenceList, int * d_zPositionReferenceList, int * d_zBackPositionReferenceList, int numberOfCoords)
{
	kdTask initialTask;
	initialTask.kdDestination = 0;
	initialTask.maxX = numberOfCoords;
	initialTask.maxY = numberOfCoords;
	initialTask.maxZ = numberOfCoords;
	initialTask.minX = 0;
	initialTask.minY = 0;
	initialTask.minZ = 0;
	taskList[0] = initialTask;
}
__global__ void online_kdConstructor_slave_nodeX(int * d_xBackPositionReferenceList, int totalNumOfConsecutiveTasks, int* kdArray, int* d_xPositionReferenceList, int*d_yPositionReferenceList, int* d_zPositionReferenceList)
{
	int id = threadIdx.x + blockDim.x * blockIdx.x;

	if (id < totalNumOfConsecutiveTasks)
	{
		kdTask task = taskList[id];
		if (task.kdDestination != -1)
		{
			int activeElementCount = 0;
			int secondCount = -1;
			int elementToBeExtracted = -1;
			int posOfExtractTarget = -1;

			for (int i = task.minX; i < task.maxX; i++) //will this loop run efficiently??
			{
				//if (d_xBackPositionReferenceList[i] != -1)
				if (d_xPositionReferenceList[d_xBackPositionReferenceList[i]] != -1) //experim
					if (d_yPositionReferenceList[d_xBackPositionReferenceList[i]] < task.maxY && d_yPositionReferenceList[d_xBackPositionReferenceList[i]] >= task.minY)
						if (d_zPositionReferenceList[d_xBackPositionReferenceList[i]] < task.maxZ && d_zPositionReferenceList[d_xBackPositionReferenceList[i]] >= task.minZ)
							activeElementCount++;
			}

			if (activeElementCount % 2 == 0)
				elementToBeExtracted = activeElementCount / 2 - 1;
			else
				elementToBeExtracted = (activeElementCount - 1) / 2;

			int i = task.minX;
			while (posOfExtractTarget == -1 && i<task.maxX)
			{
				//if (d_xBackPositionReferenceList[i] != -1)
				if (d_xPositionReferenceList[d_xBackPositionReferenceList[i]] != -1) //experimental replacement
					if (d_yPositionReferenceList[d_xBackPositionReferenceList[i]] < task.maxY && d_yPositionReferenceList[d_xBackPositionReferenceList[i]] >= task.minY)
						if (d_zPositionReferenceList[d_xBackPositionReferenceList[i]] < task.maxZ && d_zPositionReferenceList[d_xBackPositionReferenceList[i]] >= task.minZ)
							secondCount++;
				if (secondCount == elementToBeExtracted)
				{
					posOfExtractTarget = i;
					secondCount++;
				}
				i++;
			}

			kdTask lowerValTask;
			lowerValTask.kdDestination = -1;
			if (activeElementCount > 2)
			{
				lowerValTask.minX = task.minX;
				lowerValTask.maxX = posOfExtractTarget;
				lowerValTask.maxY = task.maxY;
				lowerValTask.maxZ = task.maxZ;
				lowerValTask.minY = task.minY;
				lowerValTask.minZ = task.minZ;
				lowerValTask.kdDestination = task.kdDestination * 2 + 1;
			}
			kdTask higherValTask;
			higherValTask.kdDestination = -1;
			if (activeElementCount > 1)
			{
				higherValTask.maxX = task.maxX;
				higherValTask.maxY = task.maxY;
				higherValTask.maxZ = task.maxZ;
				higherValTask.minX = posOfExtractTarget + 1;
				higherValTask.minY = task.minY;
				higherValTask.minZ = task.minZ;
				higherValTask.kdDestination = task.kdDestination * 2 + 2;
			}
			taskList[id] = lowerValTask;
			taskList[id + totalNumOfConsecutiveTasks] = higherValTask;
			if (posOfExtractTarget != -1)
			{
				kdArray[task.kdDestination] = d_xBackPositionReferenceList[posOfExtractTarget];
				d_xPositionReferenceList[d_xBackPositionReferenceList[posOfExtractTarget]] = -1;
				d_yPositionReferenceList[d_xBackPositionReferenceList[posOfExtractTarget]] = -1;
				d_zPositionReferenceList[d_xBackPositionReferenceList[posOfExtractTarget]] = -1;
			}

		}
	}
}
__global__ void online_kdConstructor_slave_nodeY(int * d_yBackPositionReferenceList, int totalNumOfConsecutiveTasks, int* kdArray, int* d_xPositionReferenceList, int*d_yPositionReferenceList, int* d_zPositionReferenceList)
{

	int id = threadIdx.x + blockDim.x * blockIdx.x;

	if (id < totalNumOfConsecutiveTasks)
	{
		kdTask task = taskList[id];
		if (task.kdDestination != -1)
		{
			int activeElementCount = 0;
			int secondCount = -1;
			int elementToBeExtracted = -1;
			int posOfExtractTarget = -1;
			for (int i = task.minY; i < task.maxY; i++) //will this loop run efficiently??
			{
				//if (d_yBackPositionReferenceList[i] != -1)
				if (d_yPositionReferenceList[d_yBackPositionReferenceList[i]] != -1)//experimental replacement
					if (d_xPositionReferenceList[d_yBackPositionReferenceList[i]] < task.maxX && d_xPositionReferenceList[d_yBackPositionReferenceList[i]] >= task.minX)
						if (d_zPositionReferenceList[d_yBackPositionReferenceList[i]] < task.maxZ && d_zPositionReferenceList[d_yBackPositionReferenceList[i]] >= task.minZ)
							activeElementCount++;
			}

			if (activeElementCount % 2 == 0)
				elementToBeExtracted = activeElementCount / 2 - 1;
			else
				elementToBeExtracted = (activeElementCount - 1) / 2;

			int i = task.minY;
			while (posOfExtractTarget == -1 && i<task.maxY)
			{
				//if (d_yBackPositionReferenceList[i] != -1)
				if (d_yPositionReferenceList[d_yBackPositionReferenceList[i]] != -1)//experimental replacement
					if (d_xPositionReferenceList[d_yBackPositionReferenceList[i]] < task.maxX && d_xPositionReferenceList[d_yBackPositionReferenceList[i]] >= task.minX)
						if (d_zPositionReferenceList[d_yBackPositionReferenceList[i]] < task.maxZ && d_zPositionReferenceList[d_yBackPositionReferenceList[i]] >= task.minZ)
							secondCount++;
				if (secondCount == elementToBeExtracted)
				{
					posOfExtractTarget = i;
					secondCount++;
				}
				i++;
			}

			kdTask lowerValTask;
			lowerValTask.kdDestination = -1;
			if (activeElementCount > 2)
			{

				lowerValTask.maxX = task.maxX;
				lowerValTask.maxY = posOfExtractTarget;
				lowerValTask.maxZ = task.maxZ;
				lowerValTask.minX = task.minX;
				lowerValTask.minY = task.minY;
				lowerValTask.minZ = task.minZ;
				lowerValTask.kdDestination = task.kdDestination * 2 + 1;
			}

			kdTask higherValTask;
			higherValTask.kdDestination = -1;
			if (activeElementCount > 1)
			{
				higherValTask.maxX = task.maxX;
				higherValTask.maxY = task.maxY;
				higherValTask.maxZ = task.maxZ;
				higherValTask.minX = task.minX;
				higherValTask.minY = posOfExtractTarget + 1;
				higherValTask.minZ = task.minZ;
				higherValTask.kdDestination = task.kdDestination * 2 + 2;
			}
			taskList[id] = lowerValTask;
			taskList[id + totalNumOfConsecutiveTasks] = higherValTask;
			if (posOfExtractTarget != -1)
			{
				kdArray[task.kdDestination] = d_yBackPositionReferenceList[posOfExtractTarget];
				d_xPositionReferenceList[d_yBackPositionReferenceList[posOfExtractTarget]] = -1;
				d_yPositionReferenceList[d_yBackPositionReferenceList[posOfExtractTarget]] = -1;
				d_zPositionReferenceList[d_yBackPositionReferenceList[posOfExtractTarget]] = -1;
			}
		}
	}
}
__global__ void online_kdConstructor_slave_nodeZ(int * d_zBackPositionReferenceList, int totalNumOfConsecutiveTasks, int* kdArray, int* d_xPositionReferenceList, int*d_yPositionReferenceList, int* d_zPositionReferenceList)
{
	int id = threadIdx.x + blockDim.x * blockIdx.x;

	if (id < totalNumOfConsecutiveTasks)
	{
		kdTask task = taskList[id];
		if (task.kdDestination != -1)
		{
			int activeElementCount = 0;
			int secondCount = -1;
			int elementToBeExtracted = -1;
			int posOfExtractTarget = -1;

			for (int i = task.minZ; i < task.maxZ; i++) //will this loop run efficiently??
			{
				//if (d_zBackPositionReferenceList[i] != -1)
				if (d_zPositionReferenceList[d_zBackPositionReferenceList[i]] != -1)//experimental replacement
					if (d_yPositionReferenceList[d_zBackPositionReferenceList[i]] < task.maxY && d_yPositionReferenceList[d_zBackPositionReferenceList[i]] >= task.minY)
						if (d_xPositionReferenceList[d_zBackPositionReferenceList[i]] < task.maxX && d_xPositionReferenceList[d_zBackPositionReferenceList[i]] >= task.minX)
							activeElementCount++;
			}

			if (activeElementCount % 2 == 0)
				elementToBeExtracted = activeElementCount / 2 - 1;
			else
				elementToBeExtracted = (activeElementCount - 1) / 2;

			int i = task.minZ;
			while (posOfExtractTarget == -1 && i<task.maxZ)
			{
				//if (d_zBackPositionReferenceList[i] != -1)
				if (d_zPositionReferenceList[d_zBackPositionReferenceList[i]] != -1)//experimental replacement
					if (d_yPositionReferenceList[d_zBackPositionReferenceList[i]] < task.maxY && d_yPositionReferenceList[d_zBackPositionReferenceList[i]] >= task.minY)
						if (d_xPositionReferenceList[d_zBackPositionReferenceList[i]] < task.maxX && d_xPositionReferenceList[d_zBackPositionReferenceList[i]] >= task.minX)
							secondCount++;
				if (secondCount == elementToBeExtracted)
				{
					posOfExtractTarget = i;
					secondCount++;
				}
				i++;
			}

			kdTask lowerValTask;
			lowerValTask.kdDestination = -1;
			if (activeElementCount > 2)
			{

				lowerValTask.maxX = task.maxX;
				lowerValTask.maxY = task.maxY;
				lowerValTask.maxZ = posOfExtractTarget;
				lowerValTask.minX = task.minX;
				lowerValTask.minY = task.minY;
				lowerValTask.minZ = task.minZ;
				lowerValTask.kdDestination = task.kdDestination * 2 + 1;
			}
			kdTask higherValTask;
			higherValTask.kdDestination = -1;
			if (activeElementCount > 1)
			{
				higherValTask.maxX = task.maxX;
				higherValTask.maxY = task.maxY;
				higherValTask.maxZ = task.maxZ;
				higherValTask.minX = task.minX;
				higherValTask.minY = task.minY;
				higherValTask.minZ = posOfExtractTarget + 1;
				higherValTask.kdDestination = task.kdDestination * 2 + 2;
			}
			taskList[id] = lowerValTask;
			taskList[id + totalNumOfConsecutiveTasks] = higherValTask;
			if (posOfExtractTarget != -1)
			{
				kdArray[task.kdDestination] = d_zBackPositionReferenceList[posOfExtractTarget];
				d_xPositionReferenceList[d_zBackPositionReferenceList[posOfExtractTarget]] = -1;
				d_yPositionReferenceList[d_zBackPositionReferenceList[posOfExtractTarget]] = -1;
				d_zPositionReferenceList[d_zBackPositionReferenceList[posOfExtractTarget]] = -1;
			};
		}
	}
}
