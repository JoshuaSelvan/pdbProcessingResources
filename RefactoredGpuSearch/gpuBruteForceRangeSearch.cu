#include"gpuBruteForceSearch.cuh"



#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
	if (code != cudaSuccess)
	{
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}



void gpuBruteForceRangeSearchAllLoadedSets(std::string atomA, std::string atomB, int requiredProximity, ProteinDataHandler heldProteinSets, AtomToNumHashTable atomReferenceTable)
{
	std::cout << "GPU BRUTE FORCE RANGE SEARCH" << std::endl;
	//all data has been preloaded into the host side ProteinDataHandler object. All that is needed is to initialise containers, moves sets of data to the gpu, process those sets of data and then return the results.

	//Initialise host and device data holders
	gpuRangeSearchResources rangeSearch;
	int safeHolderSize = 32000 * sizeof(int)*2; //home gpu probably cant hold this many -.-
	//int blocks = 0;
	//int threads = 0;
	//int currentEntry = 0;
	//gpuErrchk(cudaMalloc((void**)&rangeSearch.d_atomACount, sizeof(int)));
	//gpuErrchk(cudaMalloc((void**)&rangeSearch.d_atomAPositionList, safeHolderSize));
	//gpuErrchk(cudaMalloc((void**)&rangeSearch.d_atomACurrentSearchDimensions, safeHolderSize));
	//gpuErrchk(cudaMalloc((void**)&rangeSearch.d_atomACurrentSearchKdTreePositions, safeHolderSize));
	//gpuErrchk(cudaMalloc((void**)&rangeSearch.d_atomAMatches, safeHolderSize));
	//gpuErrchk(cudaMalloc((void**)&rangeSearch.d_atomBMatches, safeHolderSize));
	//gpuErrchk(cudaMalloc((void**)&rangeSearch.d_MatchesCount, sizeof(int)));
//	gpuErrchk(cudaMalloc((void**)&rangeSearch.d_nextSearchCount, sizeof(int)));
//	gpuErrchk(cudaMalloc((void**)&rangeSearch.d_completionFlag, sizeof(int)));
//	gpuErrchk(cudaMalloc((void**)&rangeSearch.d_kdTreeSets, safeHolderSize * 2));

	gpuErrchk(cudaMalloc((void**)&rangeSearch.d_namesSets, safeHolderSize / sizeof(int) * sizeof(short)));
	
	//gpuErrchk(cudaMalloc((void**)&rangeSearch.d_xyzCoordsSets, safeHolderSize * 3));

	gpuErrchk(cudaMalloc((void**)&rangeSearch.d_xCoordsSets, safeHolderSize));
	gpuErrchk(cudaMalloc((void**)&rangeSearch.d_yCoordsSets, safeHolderSize));
	gpuErrchk(cudaMalloc((void**)&rangeSearch.d_zCoordsSets, safeHolderSize));

//	rangeSearch.h_atomACount = (int*)malloc(sizeof(int));
//	rangeSearch.h_atomAPositionList = (int*)malloc(safeHolderSize);
//	rangeSearch.h_nextSearchCount = (int*)malloc(sizeof(int));
//	rangeSearch.h_completionFlag = (int*)malloc(sizeof(int));
//	rangeSearch.h_MatchesCount = (int*)malloc(sizeof(int));
//	rangeSearch.h_atomAMatches = (int*)malloc(safeHolderSize);
//	rangeSearch.h_atomBMatches = (int*)malloc(safeHolderSize);

	cudaStream_t streams[3];


	cudaStreamCreate(&streams[0]);
	cudaStreamCreate(&streams[1]);
	cudaStreamCreate(&streams[2]);

	//load data into holder arrays:
	int entrySize = heldProteinSets.ProteinDataHolder[0].MaxEntrySize;
	int kdEntrySize = heldProteinSets.ProteinDataHolder[0].KdTreeSize;

	//gpuErrchk(cudaMemcpyAsync(rangeSearch.d_xCoordsSets, heldProteinSets.ProteinDataHolder[0].xCoordsSets, entrySize*sizeof(int), cudaMemcpyHostToDevice, streams[0]));
	//gpuErrchk(cudaMemcpyAsync(rangeSearch.d_yCoordsSets, heldProteinSets.ProteinDataHolder[0].yCoordsSets, entrySize*sizeof(int), cudaMemcpyHostToDevice, streams[1]));
	//gpuErrchk(cudaMemcpyAsync(rangeSearch.d_zCoordsSets, heldProteinSets.ProteinDataHolder[0].zCoordsSets, entrySize*sizeof(int), cudaMemcpyHostToDevice, streams[2]));
	//Leaving this here as backup notation for the kd tree search - delete when its operational
	//gpuErrchk(cudaMemcpyAsync(rangeSearch.d_xyzCoordsSets, heldProteinSets.ProteinDataHolder[0].xCoordsSets, entrySize*sizeof(int), cudaMemcpyHostToDevice, streams[0]));
	//gpuErrchk(cudaMemcpyAsync(rangeSearch.d_xyzCoordsSets + entrySize, heldProteinSets.ProteinDataHolder[0].yCoordsSets, entrySize*sizeof(int), cudaMemcpyHostToDevice, streams[1]));
	//gpuErrchk(cudaMemcpyAsync(rangeSearch.d_xyzCoordsSets + entrySize * 2, heldProteinSets.ProteinDataHolder[0].zCoordsSets, entrySize*sizeof(int), cudaMemcpyHostToDevice, streams[2]));
	//gpuErrchk(cudaMemcpyAsync(rangeSearch.d_namesSets, heldProteinSets.ProteinDataHolder[0].namesSets, entrySize / sizeof(int)*sizeof(short), cudaMemcpyHostToDevice, streams[0]));
//	gpuErrchk(cudaMemcpyAsync(rangeSearch.d_kdTreeSets, heldProteinSets.ProteinDataHolder[0].kdTrees, kdEntrySize*sizeof(int), cudaMemcpyHostToDevice, streams[0]));


	//initiate search
	int soughtAtomANumber = 1414;
	int soughtAtomBNumber =1514;

	//int soughtAtomANumber = atomReferenceTable.retrieveHashValue(atomA);//heldProteinSets.ProteinDataHolder[0].namesSets[0];
	//int soughtAtomBNumber = atomReferenceTable.retrieveHashValue(atomB);//heldProteinSets.ProteinDataHolder[0].namesSets[1];

	//gpuErrchk(cudaMemsetAsync(rangeSearch.d_atomACount, 0, sizeof(int), streams[0]));
	//gpuErrchk(cudaMemsetAsync(rangeSearch.d_atomAPositionList, -1, safeHolderSize, streams[1]));
	//gpuErrchk(cudaMemsetAsync(rangeSearch.d_atomAMatches, -1, safeHolderSize, streams[2]));
	//gpuErrchk(cudaMemsetAsync(rangeSearch.d_atomBMatches, -1, safeHolderSize, streams[0]));
//	gpuErrchk(cudaMemsetAsync(rangeSearch.d_atomACurrentSearchDimensions, 0, safeHolderSize, streams[1]));
//	gpuErrchk(cudaMemsetAsync(rangeSearch.d_atomACurrentSearchKdTreePositions, -1, safeHolderSize, streams[2]));
//	gpuErrchk(cudaMemsetAsync(rangeSearch.d_MatchesCount, 0, sizeof(int), streams[0]));
//	gpuErrchk(cudaMemsetAsync(rangeSearch.d_nextSearchCount, 0, sizeof(int), streams[1]));
//	gpuErrchk(cudaMemsetAsync(rangeSearch.d_completionFlag, 0, sizeof(int), streams[2]));
//	rangeSearch.h_MatchesCount[0] = 0;

//	calculateInitialBlocksAndThreads(blocks, threads, heldProteinSets.ProteinDataHolder[0].MaxEntrySize);



//	cudaDeviceSynchronize();

//	device_side__locateElement << <blocks, threads, 0, streams[0] >> >(rangeSearch.d_namesSets, currentEntry, soughtAtomANumber, rangeSearch.d_atomAPositionList, rangeSearch.d_atomACount, rangeSearch.d_atomACurrentSearchDimensions, rangeSearch.d_atomACurrentSearchKdTreePositions, heldProteinSets.ProteinDataHolder[0].MaxEntrySize);

//	gpuErrchk(cudaMemcpyAsync(rangeSearch.h_atomACount, rangeSearch.d_atomACount, sizeof(int), cudaMemcpyDeviceToHost, streams[0]));
//	gpuErrchk(cudaMemcpyAsync(rangeSearch.h_atomAPositionList, rangeSearch.d_atomAPositionList, safeHolderSize, cudaMemcpyDeviceToHost, streams[1]));


	for (int i = 0; i < 5; i++)
	{
		

		std::cout << "Transfering set: " << i << " to the gpu" << std::endl;
		gpuErrchk(cudaMemcpy(rangeSearch.d_namesSets, heldProteinSets.ProteinDataHolder[i].namesSets, heldProteinSets.ProteinDataHolder[i].MaxEntrySize * heldProteinSets.ProteinDataHolder[i].heldEntries * sizeof(short), cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(rangeSearch.d_xCoordsSets, heldProteinSets.ProteinDataHolder[i].xCoordsSets, heldProteinSets.ProteinDataHolder[i].MaxEntrySize * heldProteinSets.ProteinDataHolder[i].heldEntries * sizeof(int), cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(rangeSearch.d_yCoordsSets, heldProteinSets.ProteinDataHolder[i].yCoordsSets, heldProteinSets.ProteinDataHolder[i].MaxEntrySize * heldProteinSets.ProteinDataHolder[i].heldEntries * sizeof(int), cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(rangeSearch.d_zCoordsSets, heldProteinSets.ProteinDataHolder[i].zCoordsSets, heldProteinSets.ProteinDataHolder[i].MaxEntrySize * heldProteinSets.ProteinDataHolder[i].heldEntries * sizeof(int), cudaMemcpyHostToDevice));

		//run bruteForce set runner
		bruteForceSearchPreLoadedArraySets(rangeSearch.d_namesSets, rangeSearch.d_xCoordsSets, rangeSearch.d_yCoordsSets, rangeSearch.d_zCoordsSets, heldProteinSets.ProteinDataHolder[i].MaxEntrySize, soughtAtomANumber, soughtAtomBNumber, heldProteinSets.ProteinDataHolder[i].heldEntries, requiredProximity, heldProteinSets.ProteinDataHolder[i].xCoordsSets, heldProteinSets.ProteinDataHolder[i].yCoordsSets, heldProteinSets.ProteinDataHolder[i].zCoordsSets, heldProteinSets.ProteinDataHolder[i].namesSets);
	
	}


	cudaStreamDestroy(streams[0]);
	cudaStreamDestroy(streams[1]);
	cudaStreamDestroy(streams[2]);
}




void bruteForceSearchPreLoadedArraySets(short * d_namesSet, int*d_xValsSet, int*d_yValsSet, int*d_zValsSet, int MaximumlengthOfChains, short elementA, short  elementB, int numOfEntries, int requiredProximity, int * h_xValsSets, int * h_yValSets, int *h_zValSets, short* h_names)
{
	
	gpuBruteForceResources resources;
	//declare variables required for later stages
	resources.blocks;
	resources.threads;
	resources.concurrentThreads;
	resources.h_resultsCount = (int*)malloc(1 * sizeof(int));
	resources.h_resultsA = (int*)malloc(MaximumlengthOfChains * 100 * sizeof(int));
	resources.h_resultsB = (int*)malloc(MaximumlengthOfChains * 100 * sizeof(int));
	resources.h_aCount = (int*)malloc(1 * sizeof(int));
	resources.h_bCount = (int*)malloc(1 * sizeof(int));
	resources.h_elementAList = (int*)malloc(MaximumlengthOfChains * sizeof(int)); //I used these for loading the element arrays back from the device to check what was in them, otherwise commented out.
	resources.h_elementBList = (int*)malloc(MaximumlengthOfChains * sizeof(int));
	resources.d_elementAList; //I used these for loading the element arrays back from the device to check what was in them, otherwise commented out.
	resources.d_elementBList;
	resources.d_resultsCount;
	resources.d_resultsA;
	resources.d_resultsB;
	resources.d_aCount;
	resources.d_bCount;


	gpuErrchk(cudaMalloc((void**)&resources.d_resultsCount, 1 * sizeof(int)));
	gpuErrchk(cudaMalloc((void**)&resources.d_resultsA, MaximumlengthOfChains * 100 * sizeof(int)));
	gpuErrchk(cudaMalloc((void**)&resources.d_resultsB, MaximumlengthOfChains * 100 * sizeof(int)));
	gpuErrchk(cudaMalloc((void**)&resources.d_elementAList, MaximumlengthOfChains * 10 * sizeof(int)));
	gpuErrchk(cudaMalloc((void**)&resources.d_elementBList, MaximumlengthOfChains * 10 * sizeof(int)));
	gpuErrchk(cudaMalloc((void**)&resources.d_aCount, 1 * sizeof(int)));
	gpuErrchk(cudaMalloc((void**)&resources.d_bCount, 1 * sizeof(int)));
	gpuErrchk(cudaMemset(resources.d_resultsCount, 0, 1 * sizeof(int)));
	gpuErrchk(cudaMemset(resources.d_resultsA, -1, MaximumlengthOfChains * 10 * sizeof(int)));
	gpuErrchk(cudaMemset(resources.d_resultsB, -1, MaximumlengthOfChains * 10 * sizeof(int)));
	gpuErrchk(cudaMemset(resources.d_aCount, 0, 1 * sizeof(int)));
	gpuErrchk(cudaMemset(resources.d_bCount, 0, 1 * sizeof(int)));

	//I assume these are for testing purposes
	gpuErrchk(cudaMemcpy(resources.h_resultsCount, resources.d_resultsCount, 1 * sizeof(int), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(resources.h_resultsA, resources.d_resultsA, 200 * sizeof(int), cudaMemcpyDeviceToHost));

	//Configure the block and thread size based on max chain length size.
	resources.threads = 512;
	if (MaximumlengthOfChains < 513)
	{
		resources.blocks = 1;
	}
	else
	{
		resources.blocks = (MaximumlengthOfChains + resources.threads - 1) / resources.threads; //may need to check this behaves as expected.}
	}
	resources.concurrentThreads = resources.blocks*resources.threads;


	short * names = (short*)malloc(MaximumlengthOfChains*numOfEntries*sizeof(short));

	int numberOfFilesProcessed = 0;
	for (short i = 0; i < numOfEntries; i++)
	{
		gpuErrchk(cudaMemset(resources.d_resultsCount, 0, 1 * sizeof(int))); //working
		gpuErrchk(cudaMemset(resources.d_resultsA, -1, MaximumlengthOfChains * 100 * sizeof(int)));
		gpuErrchk(cudaMemset(resources.d_resultsB, -1, MaximumlengthOfChains * 100 * sizeof(int)));
		gpuErrchk(cudaMemset(resources.d_aCount, 0, 1 * sizeof(int)));
		gpuErrchk(cudaMemset(resources.d_bCount, 0, 1 * sizeof(int)));
		gpuErrchk(cudaMemset(resources.d_elementAList, -1, MaximumlengthOfChains * 10 * sizeof(int)));
		gpuErrchk(cudaMemset(resources.d_elementBList, -1, MaximumlengthOfChains * 10 * sizeof(int)));

		resources.h_aCount[0] = 0;
		resources.h_bCount[0] = 0;




		// -- -- experimental combined kernel
	//	DeviceLoadedArrays_SingleProtein_CompleteBruteForceSearch << <resources.blocks, resources.threads >> >(d_namesSet, resources.d_elementAList, resources.d_elementBList, elementA, elementB, resources.d_aCount, resources.d_bCount, MaximumlengthOfChains, i, resources.concurrentThreads, resources.d_resultsA, resources.d_resultsB, resources.d_resultsCount, d_xValsSet, d_yValsSet, d_zValsSet, requiredProximity);
		DeviceLoadedArrays_SingleProtein_LocateElements << <resources.blocks, resources.threads >> >(d_namesSet, resources.d_elementAList, resources.d_elementBList, elementA, elementB, resources.d_aCount, resources.d_bCount, MaximumlengthOfChains, i, resources.concurrentThreads, resources.d_resultsA, resources.d_resultsB, resources.d_resultsCount, d_xValsSet, d_yValsSet, d_zValsSet, requiredProximity);
		DeviceLoadedArrays_SingleProtein_BruteForceSearch << <resources.blocks, resources.threads >> >(d_namesSet, resources.d_elementAList, resources.d_elementBList, elementA, elementB, resources.d_aCount, resources.d_bCount, MaximumlengthOfChains, i, resources.concurrentThreads, resources.d_resultsA, resources.d_resultsB, resources.d_resultsCount, d_xValsSet, d_yValsSet, d_zValsSet, requiredProximity);
		
		gpuErrchk(cudaMemcpy(resources.h_elementAList, resources.d_elementAList, MaximumlengthOfChains * sizeof(int), cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpy(resources.h_elementBList, resources.d_elementBList, MaximumlengthOfChains * sizeof(int), cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpy(resources.h_aCount, resources.d_aCount, 1 * sizeof(int), cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpy(resources.h_bCount, resources.d_bCount, 1 * sizeof(int), cudaMemcpyDeviceToHost));

		std::cout << "File: " << i << " - atom A count: " << resources.h_aCount[0] << "  ";
		//std::cout << "Number of A atoms in file " << numOfEntries << " is: " << resources.h_aCount[0] << std::endl;


		//retrieve result arrays from device
		gpuErrchk(cudaMemcpy(resources.h_resultsCount, resources.d_resultsCount, sizeof(int), cudaMemcpyDeviceToHost));

		//cout<<"Result count: "<<h_resultsCount[0]<<endl;
		if (resources.h_resultsCount[0]<(MaximumlengthOfChains * 9))
		{


			gpuErrchk(cudaMemcpy(resources.h_resultsA, resources.d_resultsA, resources.h_resultsCount[0] * 2 * sizeof(int), cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpy(resources.h_resultsB, resources.d_resultsB, resources.h_resultsCount[0] * 2 * sizeof(int), cudaMemcpyDeviceToHost));

			if (resources.h_resultsCount[0]>0)
			{

				if ((resources.h_bCount[0]>0) || (resources.h_aCount[0]>0))
				{

					int larger = resources.h_aCount[0];
					if (resources.h_bCount[0]>larger)
						larger = resources.h_bCount[0];
					for (int t = 0; t < larger; t++)
					{
						//              outputFile << h_elementAList[t]<<"\t"<<h_names[i*MaximumlengthOfChains+h_elementAList[t]] << "\t" << h_elementBList[t] << "\t"<< h_names[i*MaximumlengthOfChains+h_elementBList[t]] <<endl;
					}
				}
				if (resources.h_resultsCount[0]>0)
				{
					for (int j = 0; j < resources.h_resultsCount[0]; j++)
					{
						//std::cout << "AtomA: " << resources.h_resultsA[j] << "\t" << (double(h_xValsSets[i*MaximumlengthOfChains + resources.h_resultsA[j]]) / 1000) << "\t" << (double(h_yValSets[i*MaximumlengthOfChains + resources.h_resultsA[j]]) / 1000) << "\t" << (double(h_zValSets[i*MaximumlengthOfChains + resources.h_resultsA[j]]) / 1000) << " atomNumber: "<<  j<<std::endl;
						//std::cout << "AtomB: " << resources.h_resultsB[j] << "\t" << (double(h_xValsSets[i*MaximumlengthOfChains + resources.h_resultsB[j]]) / 1000) << "\t" << (double(h_yValSets[i*MaximumlengthOfChains + resources.h_resultsB[j]]) / 1000) << "\t" << (double(h_zValSets[i*MaximumlengthOfChains + resources.h_resultsB[j]]) / 1000) << std::endl;
						//std::cout << std::endl;
					}
				}
			}
		}
		//if (resources.h_resultsCount[0] != 0)
			std::cout << "Matches in file: " << resources.h_resultsCount[0] << std::endl;
			//std::cout << "Number of matches in file " << numberOfFilesProcessed << " is: " << resources.h_resultsCount[0] << std::endl;
		numberOfFilesProcessed++;
	}
	std::cout << "Number of entries is: " << numOfEntries << std::endl;
	gpuErrchk(cudaFree(resources.d_resultsCount));
	gpuErrchk(cudaFree(resources.d_resultsA));
	gpuErrchk(cudaFree(resources.d_resultsB));
	gpuErrchk(cudaFree(resources.d_elementAList));
	gpuErrchk(cudaFree(resources.d_elementBList));
	gpuErrchk(cudaFree(resources.d_aCount));
	gpuErrchk(cudaFree(resources.d_bCount));
};




__global__ void DeviceLoadedArrays_SingleProtein_CompleteBruteForceSearch(short * d_namesSet, int* d_elementAList, int* d_elementBList, short atomA, short atomB, int * d_aCount, int * d_bCount, int standardizedEntrySize, int entryNumber, int concurrentThreads, int*d_resultsAList, int*d_resultsBList, int*d_resultsCount, int* d_xValsSet, int* d_yValsSet, int* d_zValsSet, int requiredProximity)
{
	int id = threadIdx.x + blockDim.x * blockIdx.x;
	int secondroundId = id;
	int resultsArrayInsertPosition;
	int currentDeviceName;
	int currentOffset = entryNumber*standardizedEntrySize;


	if (id == 0)
	{
		d_resultsCount[0] = 0;
		d_aCount[0] = 0;
		d_bCount[0] = 0;
	}
	__syncthreads();
	//while (id < (standardizedEntrySize +currentOffset))
	{
		currentDeviceName = d_namesSet[id + currentOffset];
		if (currentDeviceName == atomA)
		{
			resultsArrayInsertPosition = atomicAdd(&d_aCount[0], 1);
			d_elementAList[resultsArrayInsertPosition] = id;
			
		}
		if (currentDeviceName == atomB)
		{
			resultsArrayInsertPosition = atomicAdd(&d_bCount[0], 1);
			d_elementBList[resultsArrayInsertPosition] = id;
		}
	//	id = id + concurrentThreads;
	}
	__syncthreads();



	if (secondroundId<d_aCount[0])
	{
		int requiredProximitySquared = requiredProximity*requiredProximity;

		int insertPosition;
		short localAtomA = d_elementAList[secondroundId];

		int currentSetOffset = entryNumber*standardizedEntrySize;

		int localXCoord = d_xValsSet[currentSetOffset + localAtomA];
		int localYCoord = d_yValsSet[currentSetOffset + localAtomA];
		int localZCoord = d_zValsSet[currentSetOffset + localAtomA];

		short currentAtomB;
		int BCurrentXCoord;
		int BCurrentYCoord;
		int BCurrentZCoord;
		int distanceBetweenAtoms;

		for (int i = 0; i < d_bCount[0]; i++)
		{
			currentAtomB = d_elementBList[i];
			BCurrentXCoord = d_xValsSet[currentSetOffset + currentAtomB];
			BCurrentYCoord = d_yValsSet[currentSetOffset + currentAtomB];
			BCurrentZCoord = d_zValsSet[currentSetOffset + currentAtomB];		
			distanceBetweenAtoms = ((localXCoord - BCurrentXCoord)*(localXCoord - BCurrentXCoord) + (localYCoord - BCurrentYCoord)*(localYCoord - BCurrentYCoord) + (localZCoord - BCurrentZCoord)*(localZCoord - BCurrentZCoord));
			
			int f = (requiredProximitySquared-distanceBetweenAtoms)/100;
		
			//This loop is triggering regardless of it being right or wrong -.-
			//if (distanceBetweenAtoms < requiredProximitySquared)
			if (f>0)
			{
				insertPosition = atomicAdd(&d_resultsCount[0], 1);
				d_resultsAList[insertPosition] = localAtomA;
				d_resultsBList[insertPosition] = currentAtomB;
			}
		}
	}
};





__global__ void DeviceLoadedArrays_SingleProtein_LocateElements(short * d_namesSet, int* d_elementAList, int* d_elementBList, short atomA, short atomB, int * d_aCount, int * d_bCount, int standardizedEntrySize, int entryNumber, int concurrentThreads, int*d_resultsAList, int*d_resultsBList, int*d_resultsCount, int* d_xValsSet, int* d_yValsSet, int* d_zValsSet, int requiredProximity)
{
	int id = threadIdx.x + blockDim.x * blockIdx.x;
	int secondroundId = id;
	int resultsArrayInsertPosition;
	int currentDeviceName;
	int currentOffset = entryNumber*standardizedEntrySize;


	if (id == 0)
	{
		d_resultsCount[0] = 0;
		d_aCount[0] = 0;
		d_bCount[0] = 0;
	}
	__syncthreads();
	//while (id < (standardizedEntrySize +currentOffset))
	{
		currentDeviceName = d_namesSet[id + currentOffset];
		if (currentDeviceName == atomA)
		{
			resultsArrayInsertPosition = atomicAdd(&d_aCount[0], 1);
			d_elementAList[resultsArrayInsertPosition] = id;

		}
		if (currentDeviceName == atomB)
		{
			resultsArrayInsertPosition = atomicAdd(&d_bCount[0], 1);
			d_elementBList[resultsArrayInsertPosition] = id;
		}
		//	id = id + concurrentThreads;
	}

};



__global__ void DeviceLoadedArrays_SingleProtein_BruteForceSearch(short * d_namesSet, int* d_elementAList, int* d_elementBList, short atomA, short atomB, int * d_aCount, int * d_bCount, int standardizedEntrySize, int entryNumber, int concurrentThreads, int*d_resultsAList, int*d_resultsBList, int*d_resultsCount, int* d_xValsSet, int* d_yValsSet, int* d_zValsSet, int requiredProximity)
{
	int secondroundId = threadIdx.x + blockDim.x * blockIdx.x;
	

	if (secondroundId<d_aCount[0])
	{
		int requiredProximitySquared = requiredProximity*requiredProximity;

		int insertPosition;
		short localAtomA = d_elementAList[secondroundId];

		int currentSetOffset = entryNumber*standardizedEntrySize;

		int localXCoord = d_xValsSet[currentSetOffset + localAtomA];
		int localYCoord = d_yValsSet[currentSetOffset + localAtomA];
		int localZCoord = d_zValsSet[currentSetOffset + localAtomA];

		short currentAtomB;
		int BCurrentXCoord;
		int BCurrentYCoord;
		int BCurrentZCoord;
		int distanceBetweenAtoms;

		for (int i = 0; i < d_bCount[0]; i++)
		{
			currentAtomB = d_elementBList[i];
			BCurrentXCoord = d_xValsSet[currentSetOffset + currentAtomB];
			BCurrentYCoord = d_yValsSet[currentSetOffset + currentAtomB];
			BCurrentZCoord = d_zValsSet[currentSetOffset + currentAtomB];
			distanceBetweenAtoms = ((localXCoord - BCurrentXCoord)*(localXCoord - BCurrentXCoord) + (localYCoord - BCurrentYCoord)*(localYCoord - BCurrentYCoord) + (localZCoord - BCurrentZCoord)*(localZCoord - BCurrentZCoord));

			int f = (requiredProximitySquared - distanceBetweenAtoms) / 100;

			//This loop is triggering regardless of it being right or wrong -.-
			//if (distanceBetweenAtoms < requiredProximitySquared)
			if (f>0)
			{
				insertPosition = atomicAdd(&d_resultsCount[0], 1);
				d_resultsAList[insertPosition] = localAtomA;
				d_resultsBList[insertPosition] = currentAtomB;
			}
		}
	}
};
