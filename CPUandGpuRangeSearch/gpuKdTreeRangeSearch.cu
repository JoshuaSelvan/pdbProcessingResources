#include"gpuKdTreeRangeSearch.cuh"

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
	if (code != cudaSuccess)
	{
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}


__global__ void doSillyAdditionsOnTheGpu(int* d_array)
{
	for (int i = 0; i < 7; i++)
	{

		d_array[i] = d_array[i] + 5;


	}



};

__global__ void device_side__locateElement(short * d_names, int entryNo, int soughtAtomNum, int* d_soughtAtomPositionList, int* d_soughtAtomCount, int*d_kdSearchCurrentDimensionList, int*d_kdSearchCurrentTreePos, int lengthOfChain)
{
	int insertPosition;
	int atomA = soughtAtomNum;

	int id = threadIdx.x + blockDim.x * blockIdx.x;

	__syncthreads();
	if (id < lengthOfChain)
	{
		if (d_names[id + entryNo*lengthOfChain] == atomA)
		{
			//int tempdisplayer = id;
			insertPosition = atomicAdd(&d_soughtAtomCount[0], 1); //apparently atomicAdd returns the value of the variable being increased before it is increased.
			//int tempdisplayer2 = insertPosition;
			d_soughtAtomPositionList[insertPosition] = id;
			d_kdSearchCurrentDimensionList[insertPosition] = 0;
			d_kdSearchCurrentTreePos[insertPosition] = 0;
		}
	}

	return;
}

__global__ void viewIntArrayOnGPU(int* d_array)
{
	int x = 4;
	for (int i = 0; i < 7; i++)
	{
		x = x + 3;
	}
};

__global__ void viewShortArrayOnGPU(short* d_array)
{
	int x = 4;
	for (int i = 0; i < 7; i++)
	{
		x = x + 3;
	}
};


__global__ void viewCPUIntOnGPU(int anInt)
{
	int x = 4;
	for (int i = 0; i < 7; i++)
	{
		x = x + 3;
	}
};

void kdTreeSearchAllAvailableFilesSeveralGpuLoadedEntrysAtATime(rangeSearchSettings& settings, ProteinDataHandler heldProteinSets, AtomToNumHashTable atomReferenceTable)
{
	int numberOfFilesInBatchPerRange[5];// = { 2, 2, 2, 2, 1 };
	loadMultiBatchRangeSizes("MultiBatchSizes.txt", settings);

	int requiredSizeForBatchesPerRange[5];
	for (int i = 0; i < 5; i++)
	{
		numberOfFilesInBatchPerRange[i] = settings.multiBatchRangeSizes[i];
		requiredSizeForBatchesPerRange[i] = numberOfFilesInBatchPerRange[i] * heldProteinSets.ProteinDataHolder[i].MaxEntrySize;
	}

	int largestRequiredBatchSizePos = largestSizeInArrayPos(requiredSizeForBatchesPerRange, 5);
	int largestCpuBatchSize = requiredSizeForBatchesPerRange[largestRequiredBatchSizePos];

	//all data has been preloaded into the host side ProteinDataHandler object. All that is needed is to initialise containers, moves sets of data to the gpu, process those sets of data and then return the results.
	std::cout << "PERFORMING TYPE 8 RANGE SEARCH: SMALL BATCH LOAD GPU KD TREE" << std::endl;

	cudaStream_t streams[2];
	cudaStreamCreate(&streams[0]);
	cudaStreamCreate(&streams[1]);

	//Initialise host and device data holders
	int SetBeingProcessed = 0;
	int SetBeingLoaded = 1;
	gpuRangeSearchResources rangeSearchSlots[2];
	//gpuBruteForceSingleEntryResources rangeSearchSlotB;



	int IndividualEntryHolderSize = 16390 * 6; //The memory required to hold a max size coordinate array + 6 more atoms


	int blocks = 1;
	int threads = 1;
//	int currentEntry = 0;
	for (int i = 0; i < 2; i++){
		gpuErrchk(cudaMalloc((void**)&rangeSearchSlots[i].d_atomACount, sizeof(int)));
		gpuErrchk(cudaMalloc((void**)&rangeSearchSlots[i].d_atomAPositionList, IndividualEntryHolderSize * 30));
		gpuErrchk(cudaMalloc((void**)&rangeSearchSlots[i].d_atomACurrentSearchDimensions, IndividualEntryHolderSize * 30));
		gpuErrchk(cudaMalloc((void**)&rangeSearchSlots[i].d_atomACurrentSearchKdTreePositions, IndividualEntryHolderSize * 30));
		gpuErrchk(cudaMalloc((void**)&rangeSearchSlots[i].d_atomAMatches, IndividualEntryHolderSize * 5*20));
		gpuErrchk(cudaMalloc((void**)&rangeSearchSlots[i].d_atomBMatches, IndividualEntryHolderSize * 5*20));
		gpuErrchk(cudaMalloc((void**)&rangeSearchSlots[i].d_MatchesCount, sizeof(int)));
		gpuErrchk(cudaMalloc((void**)&rangeSearchSlots[i].d_nextSearchCount, sizeof(int)));
		gpuErrchk(cudaMalloc((void**)&rangeSearchSlots[i].d_completionFlag, sizeof(int)));
		gpuErrchk(cudaMalloc((void**)&rangeSearchSlots[i].d_kdTreeSets, largestCpuBatchSize*sizeof(int) * 2));
		gpuErrchk(cudaMalloc((void**)&rangeSearchSlots[i].d_namesSets, largestCpuBatchSize*sizeof(int) / 2));
		gpuErrchk(cudaMalloc((void**)&rangeSearchSlots[i].d_xyzCoordsSets, largestCpuBatchSize*sizeof(int) * 3));
		rangeSearchSlots[i].h_atomACount = (int*)malloc(sizeof(int));
		rangeSearchSlots[i].h_atomAPositionList = (int*)malloc(IndividualEntryHolderSize * 30);
		rangeSearchSlots[i].h_nextSearchCount = (int*)malloc(sizeof(int));
		rangeSearchSlots[i].h_completionFlag = (int*)malloc(sizeof(int));
		rangeSearchSlots[i].h_MatchesCount = (int*)malloc(sizeof(int));
		rangeSearchSlots[i].h_atomAMatches = (int*)malloc(IndividualEntryHolderSize *5* 30);
		rangeSearchSlots[i].h_atomBMatches = (int*)malloc(IndividualEntryHolderSize *5* 30);
	}


	int soughtAtomANumber = atomReferenceTable.retrieveHashValue(settings.AtomTypeOne);
	int soughtAtomBNumber = atomReferenceTable.retrieveHashValue(settings.AtomTypeTwo);

	int currentMaxEntrySize;
	int currentHeldEntries;

	for (int i = 0; i < 5; i++) //For each of the 5 range lengths of stored protein: Pin memory.
	{
		currentMaxEntrySize = heldProteinSets.ProteinDataHolder[i].MaxEntrySize;
		currentHeldEntries = heldProteinSets.ProteinDataHolder[i].heldEntries;
		if (currentHeldEntries > 0)
		{
			gpuErrchk(cudaHostRegister(heldProteinSets.ProteinDataHolder[i].kdTrees, currentMaxEntrySize * currentHeldEntries * sizeof(int) * 2, 0));
			gpuErrchk(cudaHostRegister(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets, currentMaxEntrySize * currentHeldEntries * sizeof(int) * 3, 0));
			gpuErrchk(cudaHostRegister(heldProteinSets.ProteinDataHolder[i].namesSets, currentMaxEntrySize * currentHeldEntries * sizeof(short), 0));
		}
	}




	for (int i = 0; i < 5; i++)//For each of the 5 range lengths of stored protein:
	{
		currentMaxEntrySize = heldProteinSets.ProteinDataHolder[i].MaxEntrySize;
		int TotalEntries = heldProteinSets.ProteinDataHolder[i].heldEntries;
		
		                int sizeOfCurrentBatches;
                if (requiredSizeForBatchesPerRange[i]>=heldProteinSets.ProteinDataHolder[i].heldEntries)
                        sizeOfCurrentBatches=requiredSizeForBatchesPerRange[i];
                else
                        sizeOfCurrentBatches=heldProteinSets.ProteinDataHolder[i].heldEntries;


		

		//int sizeOfCurrentBatches = requiredSizeForBatchesPerRange[i];
		int NumberOfBatchesToProcess = TotalEntries / numberOfFilesInBatchPerRange[i]; //Needs to be rounded up - must check
		
		clock_t n, m;


		int outputType = settings.resultsPrintFormat;
		outputHandler filePrinter;
		std::string printType;
		if (outputType == 3) { printType = "_Summary"; }
		else if (outputType == 4) { printType = "_Detailed"; }

		if (outputType == 3 || outputType == 4)	{ filePrinter.initializeOutputfile("GpuKdTreeResults_Range_", currentMaxEntrySize, "_Files_", TotalEntries, printType); }

		n = clock();

		if (TotalEntries > 0)
		{
			std::cout << "Processing Range set: " << i << std::endl;
			std::cout << "Number of present entries is: " << heldProteinSets.ProteinDataHolder[i].heldEntries << std::endl;
			for (int currentBatch = 0; currentBatch < NumberOfBatchesToProcess + 1; currentBatch++)
			{
				if (currentBatch == 0)
				{

					//Load first set of details onto the gpu but do not process them -- needs work
					gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingLoaded].d_xyzCoordsSets, heldProteinSets.ProteinDataHolder[i].xyzCoordsSets + sizeOfCurrentBatches*currentBatch * 3, sizeOfCurrentBatches*sizeof(int) * 3, cudaMemcpyHostToDevice, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingLoaded].d_namesSets, heldProteinSets.ProteinDataHolder[i].namesSets + sizeOfCurrentBatches*currentBatch, sizeOfCurrentBatches *sizeof(short), cudaMemcpyHostToDevice, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingLoaded].d_kdTreeSets, heldProteinSets.ProteinDataHolder[i].kdTrees + sizeOfCurrentBatches*currentBatch * 2, sizeOfCurrentBatches * 2 * sizeof(int), cudaMemcpyHostToDevice, streams[SetBeingLoaded]));

					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_atomACount, 0, sizeof(int), streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_atomAPositionList, -1, IndividualEntryHolderSize * 5, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_atomAMatches, -1, IndividualEntryHolderSize * 5, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_atomBMatches, -1, IndividualEntryHolderSize * 5, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_atomACurrentSearchDimensions, 0, IndividualEntryHolderSize * 5, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_atomACurrentSearchKdTreePositions, -1, IndividualEntryHolderSize * 5, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_MatchesCount, 0, sizeof(int), streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_nextSearchCount, 0, sizeof(int), streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_completionFlag, 0, sizeof(int), streams[SetBeingLoaded]));
					rangeSearchSlots[SetBeingLoaded].h_MatchesCount[0] = 0;
				}
				else if (currentBatch == NumberOfBatchesToProcess)
				{
					//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					//process all proteins in loaded batch:

					for (int currentEntry = 0; currentEntry < numberOfFilesInBatchPerRange[i]; currentEntry++)
					{
						if(currentEntry%100==0)
							std::cout<<"i";
						//current loaded protein
						gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingProcessed].d_atomACount, 0, sizeof(int), streams[SetBeingProcessed]));
						gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingProcessed].d_atomAPositionList, -1, IndividualEntryHolderSize * 5, streams[SetBeingProcessed]));
						gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingProcessed].d_atomAMatches, -1, IndividualEntryHolderSize * 5, streams[SetBeingProcessed]));
						gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingProcessed].d_atomBMatches, -1, IndividualEntryHolderSize * 5, streams[SetBeingProcessed]));
						gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingProcessed].d_atomACurrentSearchDimensions, 0, IndividualEntryHolderSize * 5, streams[SetBeingProcessed]));
						gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingProcessed].d_atomACurrentSearchKdTreePositions, -1, IndividualEntryHolderSize * 5, streams[SetBeingProcessed]));
						gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingProcessed].d_MatchesCount, 0, sizeof(int), streams[SetBeingProcessed]));
						gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingProcessed].d_nextSearchCount, 0, sizeof(int), streams[SetBeingProcessed]));
						gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingProcessed].d_completionFlag, 0, sizeof(int), streams[SetBeingProcessed]));
						rangeSearchSlots[SetBeingLoaded].h_MatchesCount[0] = 0;



						threads = 512;
						blocks = heldProteinSets.ProteinDataHolder[i].MaxEntrySize / 512;
						device_side__locateElement << <blocks, threads, 0, streams[SetBeingProcessed] >> >(rangeSearchSlots[SetBeingProcessed].d_namesSets, currentEntry, soughtAtomANumber, rangeSearchSlots[SetBeingProcessed].d_atomAPositionList, rangeSearchSlots[SetBeingProcessed].d_atomACount, rangeSearchSlots[SetBeingProcessed].d_atomACurrentSearchDimensions, rangeSearchSlots[SetBeingProcessed].d_atomACurrentSearchKdTreePositions, heldProteinSets.ProteinDataHolder[i].MaxEntrySize);

						gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingProcessed].h_atomACount, rangeSearchSlots[SetBeingProcessed].d_atomACount, sizeof(int), cudaMemcpyDeviceToHost, streams[SetBeingProcessed]));

						int maxNumberoFCycles = 13 + (heldProteinSets.ProteinDataHolder[i].MaxEntrySize) / 1024;
						rangeSearchSlots[SetBeingProcessed].blocks = 1;
						rangeSearchSlots[SetBeingProcessed].threads = 1;
						if (rangeSearchSlots[SetBeingProcessed].h_atomACount[0] > 0)
						{
							for (int p = 0; p < maxNumberoFCycles; p++)
							{

								gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingProcessed].h_atomACount, rangeSearchSlots[SetBeingProcessed].d_atomACount, sizeof(int), cudaMemcpyDeviceToHost, streams[SetBeingProcessed]));

								if (rangeSearchSlots[SetBeingProcessed].h_atomACount[0] < 513)
								{
									rangeSearchSlots[SetBeingProcessed].blocks = 1;
									rangeSearchSlots[SetBeingProcessed].threads = rangeSearchSlots[SetBeingProcessed].h_atomACount[0];
								}
								else
								{
									rangeSearchSlots[SetBeingProcessed].threads = 512;
									rangeSearchSlots[SetBeingProcessed].blocks = rangeSearchSlots[SetBeingProcessed].h_atomACount[0] / 512 + 1;
								}



								device_side_ProcessCurrentTreePositionsV3 << <rangeSearchSlots[SetBeingProcessed].blocks, rangeSearchSlots[SetBeingProcessed].threads, 0, streams[SetBeingProcessed] >> >(rangeSearchSlots[SetBeingProcessed].d_xyzCoordsSets, rangeSearchSlots[SetBeingProcessed].d_namesSets, rangeSearchSlots[SetBeingProcessed].d_kdTreeSets, rangeSearchSlots[SetBeingProcessed].d_atomAMatches, rangeSearchSlots[SetBeingProcessed].d_atomBMatches, rangeSearchSlots[SetBeingProcessed].d_MatchesCount, rangeSearchSlots[SetBeingProcessed].d_atomACurrentSearchKdTreePositions, rangeSearchSlots[SetBeingProcessed].d_atomAPositionList, rangeSearchSlots[SetBeingProcessed].d_atomACurrentSearchDimensions, rangeSearchSlots[SetBeingProcessed].d_atomACount, heldProteinSets.ProteinDataHolder[i].MaxEntrySize, settings.requiredProximity*settings.requiredProximity, soughtAtomBNumber, rangeSearchSlots[SetBeingProcessed].d_nextSearchCount, currentEntry, heldProteinSets.ProteinDataHolder[i].KdTreeSize, settings.requiredProximity, p);
								SetNextCountAsCurrentAndCheckFlag << <1, 1, 0, streams[SetBeingProcessed] >> >(rangeSearchSlots[SetBeingProcessed].d_atomACount, rangeSearchSlots[SetBeingProcessed].d_nextSearchCount, rangeSearchSlots[SetBeingProcessed].d_completionFlag); //update device side counters reflecting how many search items exist in the list


								gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingProcessed].h_completionFlag, rangeSearchSlots[SetBeingProcessed].d_completionFlag, sizeof(int), cudaMemcpyDeviceToHost, streams[SetBeingProcessed]));
								if (rangeSearchSlots[SetBeingProcessed].h_completionFlag[0] == 1)
								{
									p = 1000;
								}

							}
						}

						gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingProcessed].h_MatchesCount, rangeSearchSlots[SetBeingProcessed].d_MatchesCount, sizeof(int), cudaMemcpyDeviceToHost, streams[SetBeingProcessed]));


						if (outputType == 1 || outputType == 2)
							std::cout << "Number of matches in file " << (currentEntry + (currentBatch - 1)*numberOfFilesInBatchPerRange[i]) << " in set " << i << " is: " << rangeSearchSlots[SetBeingProcessed].h_MatchesCount[0] << std::endl;
						else if (outputType == 3 || outputType == 4)
							filePrinter.printLineToOutputFile("Number of matches in file ", (currentEntry + (currentBatch - 1)*numberOfFilesInBatchPerRange[i]), " in set ", i, "  is: ", rangeSearchSlots[SetBeingProcessed].h_MatchesCount[0]);





						if ((rangeSearchSlots[SetBeingProcessed].h_MatchesCount[0] > 0) && (rangeSearchSlots[SetBeingProcessed].h_atomACount[0] > 0))
						{
							gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingProcessed].h_atomAMatches, rangeSearchSlots[SetBeingProcessed].d_atomAMatches, sizeof(int)*rangeSearchSlots[SetBeingProcessed].h_MatchesCount[0], cudaMemcpyDeviceToHost, streams[SetBeingProcessed]));
							gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingProcessed].h_atomBMatches, rangeSearchSlots[SetBeingProcessed].d_atomBMatches, sizeof(int)*rangeSearchSlots[SetBeingProcessed].h_MatchesCount[0], cudaMemcpyDeviceToHost, streams[SetBeingProcessed]));
							int valuesOffset = heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 3 * currentEntry;

							for (int j = 0; j < rangeSearchSlots[SetBeingProcessed].h_MatchesCount[0]; j++)
							{
								if (outputType == 2)
								{
									std::cout << "Atom A: " << heldProteinSets.ProteinDataHolder[i].namesSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]] << "\t - Pos : " << rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j] << "\t X: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]])) / 1000) << "\t Y: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]])) / 1000) << "\t Z: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 2 + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]])) / 1000) << std::endl;
									std::cout << "Atom B: " << heldProteinSets.ProteinDataHolder[i].namesSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]] << "\t - Pos : " << rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j] << "\t X: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]])) / 1000) << "\t Y: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]])) / 1000) << "\t Z: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 2 + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]])) / 1000) << std::endl << std::endl;

									std::cout << std::endl;
								}
								else if (outputType == 4)
								{
									filePrinter.printLineToOutputFile("Atom A: ", heldProteinSets.ProteinDataHolder[i].namesSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]], "\t - Pos : ", rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j], "\t X: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]])) / 1000), "\t Y: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]])) / 1000), "\t Z: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 2 + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]])) / 1000));
									filePrinter.printLineToOutputFile("Atom B: ", heldProteinSets.ProteinDataHolder[i].namesSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]], "\t - Pos : ", rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j], "\t X: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]])) / 1000), "\t Y: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]])) / 1000), "\t Z: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 2 + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]])) / 1000));
									filePrinter.printLineToOutputFile("");
								}

							}
						}
						//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


					}
				}
				else
				{


					//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					//Load next batch of entries onto the gpu.

					gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingLoaded].d_xyzCoordsSets, heldProteinSets.ProteinDataHolder[i].xyzCoordsSets + sizeOfCurrentBatches*currentBatch * 3, sizeOfCurrentBatches*sizeof(int) * 3, cudaMemcpyHostToDevice, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingLoaded].d_namesSets, heldProteinSets.ProteinDataHolder[i].namesSets + sizeOfCurrentBatches*currentBatch, sizeOfCurrentBatches *sizeof(short), cudaMemcpyHostToDevice, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingLoaded].d_kdTreeSets, heldProteinSets.ProteinDataHolder[i].kdTrees + sizeOfCurrentBatches*currentBatch * 2, sizeOfCurrentBatches * 2 * sizeof(int), cudaMemcpyHostToDevice, streams[SetBeingLoaded]));

					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_atomACount, 0, sizeof(int), streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_atomAPositionList, -1, IndividualEntryHolderSize * 5, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_atomAMatches, -1, IndividualEntryHolderSize * 5, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_atomBMatches, -1, IndividualEntryHolderSize * 5, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_atomACurrentSearchDimensions, 0, IndividualEntryHolderSize * 5, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_atomACurrentSearchKdTreePositions, -1, IndividualEntryHolderSize * 5, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_MatchesCount, 0, sizeof(int), streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_nextSearchCount, 0, sizeof(int), streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_completionFlag, 0, sizeof(int), streams[SetBeingLoaded]));
					rangeSearchSlots[SetBeingLoaded].h_MatchesCount[0] = 0;


					//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					//process all proteins in loaded batch:


					for (int currentEntry = 0; currentEntry < numberOfFilesInBatchPerRange[i]; currentEntry++)
					{
						//std::cout<<"Processing entry: "<<currentEntry<<std::endl;
						//current loaded protein
						gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingProcessed].d_atomACount, 0, sizeof(int), streams[SetBeingProcessed]));
						gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingProcessed].d_atomAPositionList, -1, IndividualEntryHolderSize * 5, streams[SetBeingProcessed]));
						gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingProcessed].d_atomAMatches, -1, IndividualEntryHolderSize * 5, streams[SetBeingProcessed]));
						gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingProcessed].d_atomBMatches, -1, IndividualEntryHolderSize * 5, streams[SetBeingProcessed]));
						gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingProcessed].d_atomACurrentSearchDimensions, 0, IndividualEntryHolderSize * 5, streams[SetBeingProcessed]));
						gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingProcessed].d_atomACurrentSearchKdTreePositions, -1, IndividualEntryHolderSize * 5, streams[SetBeingProcessed]));
						gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingProcessed].d_MatchesCount, 0, sizeof(int), streams[SetBeingProcessed]));
						gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingProcessed].d_nextSearchCount, 0, sizeof(int), streams[SetBeingProcessed]));
						gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingProcessed].d_completionFlag, 0, sizeof(int), streams[SetBeingProcessed]));
						rangeSearchSlots[SetBeingLoaded].h_MatchesCount[0] = 0;



						threads = 512;
						blocks = heldProteinSets.ProteinDataHolder[i].MaxEntrySize / 512;
						device_side__locateElement << <blocks, threads, 0, streams[SetBeingProcessed] >> >(rangeSearchSlots[SetBeingProcessed].d_namesSets, currentEntry, soughtAtomANumber, rangeSearchSlots[SetBeingProcessed].d_atomAPositionList, rangeSearchSlots[SetBeingProcessed].d_atomACount, rangeSearchSlots[SetBeingProcessed].d_atomACurrentSearchDimensions, rangeSearchSlots[SetBeingProcessed].d_atomACurrentSearchKdTreePositions, heldProteinSets.ProteinDataHolder[i].MaxEntrySize);

						gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingProcessed].h_atomACount, rangeSearchSlots[SetBeingProcessed].d_atomACount, sizeof(int), cudaMemcpyDeviceToHost, streams[SetBeingProcessed]));

						int maxNumberoFCycles = 13 + (heldProteinSets.ProteinDataHolder[i].MaxEntrySize) / 1024;
						rangeSearchSlots[SetBeingProcessed].blocks = 1;
						rangeSearchSlots[SetBeingProcessed].threads = 1;
						if (rangeSearchSlots[SetBeingProcessed].h_atomACount[0] > 0)
						{
							for (int p = 0; p < maxNumberoFCycles; p++)
							{

								gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingProcessed].h_atomACount, rangeSearchSlots[SetBeingProcessed].d_atomACount, sizeof(int), cudaMemcpyDeviceToHost, streams[SetBeingProcessed]));

								if (rangeSearchSlots[SetBeingProcessed].h_atomACount[0] < 513)
								{
									rangeSearchSlots[SetBeingProcessed].blocks = 1;
									rangeSearchSlots[SetBeingProcessed].threads = rangeSearchSlots[SetBeingProcessed].h_atomACount[0];
								}
								else
								{
									rangeSearchSlots[SetBeingProcessed].threads = 512;
									rangeSearchSlots[SetBeingProcessed].blocks = rangeSearchSlots[SetBeingProcessed].h_atomACount[0] / 512 + 1;
								}


								//std::cout << "Details before print: " << std::endl;
								//std::cout << "Threads: " << rangeSearchSlots[SetBeingProcessed].threads<< std::endl;
								//std::cout << "Blocks: " << rangeSearchSlots[SetBeingProcessed].blocks<<std::endl;
								//std::cout << "total processing: " << rangeSearchSlots[SetBeingProcessed].threads*rangeSearchSlots[SetBeingProcessed].blocks << std::endl;
								//std::cout << "Matches holder size: " << IndividualEntryHolderSize * 5*10 << std::endl;
								//std::cout << "Type 1 atoms: " << rangeSearchSlots[SetBeingProcessed].h_atomACount[0]<< std::endl;
								//std::cout << "Available search slots: "<<16390<<std::endl;
								//std::cout << "active search slots:"<<

								device_side_ProcessCurrentTreePositionsV3 << <rangeSearchSlots[SetBeingProcessed].blocks, rangeSearchSlots[SetBeingProcessed].threads, 0, streams[SetBeingProcessed] >> >(rangeSearchSlots[SetBeingProcessed].d_xyzCoordsSets, rangeSearchSlots[SetBeingProcessed].d_namesSets, rangeSearchSlots[SetBeingProcessed].d_kdTreeSets, rangeSearchSlots[SetBeingProcessed].d_atomAMatches, rangeSearchSlots[SetBeingProcessed].d_atomBMatches, rangeSearchSlots[SetBeingProcessed].d_MatchesCount, rangeSearchSlots[SetBeingProcessed].d_atomACurrentSearchKdTreePositions, rangeSearchSlots[SetBeingProcessed].d_atomAPositionList, rangeSearchSlots[SetBeingProcessed].d_atomACurrentSearchDimensions, rangeSearchSlots[SetBeingProcessed].d_atomACount, heldProteinSets.ProteinDataHolder[i].MaxEntrySize, settings.requiredProximity*settings.requiredProximity, soughtAtomBNumber, rangeSearchSlots[SetBeingProcessed].d_nextSearchCount, currentEntry, heldProteinSets.ProteinDataHolder[i].KdTreeSize, settings.requiredProximity, p);

								//								gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingProcessed].h_MatchesCount, rangeSearchSlots[SetBeingProcessed].d_MatchesCount, sizeof(int), cudaMemcpyDeviceToHost, streams[SetBeingProcessed]));
								//std::cout << "Results count: " << rangeSearchSlots[SetBeingProcessed].h_MatchesCount[0]<<std::endl;
								//std::cout<<std::endl<<std::endl;
								
								SetNextCountAsCurrentAndCheckFlag << <1, 1, 0, streams[SetBeingProcessed] >> >(rangeSearchSlots[SetBeingProcessed].d_atomACount, rangeSearchSlots[SetBeingProcessed].d_nextSearchCount, rangeSearchSlots[SetBeingProcessed].d_completionFlag); //update device side counters reflecting how many search items exist in the list


								gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingProcessed].h_completionFlag, rangeSearchSlots[SetBeingProcessed].d_completionFlag, sizeof(int), cudaMemcpyDeviceToHost, streams[SetBeingProcessed]));
								if (rangeSearchSlots[SetBeingProcessed].h_completionFlag[0] == 1)
								{
									p = 1000;
								}

							}
						}

						gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingProcessed].h_MatchesCount, rangeSearchSlots[SetBeingProcessed].d_MatchesCount, sizeof(int), cudaMemcpyDeviceToHost, streams[SetBeingProcessed]));


						if (outputType == 1 || outputType == 2)
							std::cout << "Number of matches in file " << (currentEntry + (currentBatch - 1)*numberOfFilesInBatchPerRange[i]) << " in set " << i << " is: " << rangeSearchSlots[SetBeingProcessed].h_MatchesCount[0] << std::endl;
						else if (outputType == 3 || outputType == 4)
							filePrinter.printLineToOutputFile("Number of matches in file ", (currentEntry + (currentBatch - 1)*numberOfFilesInBatchPerRange[i]), " in set ", i, "  is: ", rangeSearchSlots[SetBeingProcessed].h_MatchesCount[0]);





						if ((rangeSearchSlots[SetBeingProcessed].h_MatchesCount[0] > 0) && (rangeSearchSlots[SetBeingProcessed].h_atomACount[0] > 0))
						{
							gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingProcessed].h_atomAMatches, rangeSearchSlots[SetBeingProcessed].d_atomAMatches, sizeof(int)*rangeSearchSlots[SetBeingProcessed].h_MatchesCount[0], cudaMemcpyDeviceToHost, streams[SetBeingProcessed]));
							gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingProcessed].h_atomBMatches, rangeSearchSlots[SetBeingProcessed].d_atomBMatches, sizeof(int)*rangeSearchSlots[SetBeingProcessed].h_MatchesCount[0], cudaMemcpyDeviceToHost, streams[SetBeingProcessed]));
							int valuesOffset = heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 3 * currentEntry;

							for (int j = 0; j < rangeSearchSlots[SetBeingProcessed].h_MatchesCount[0]; j++)
							{
								if (outputType == 2)
								{
									std::cout << "Atom A: " << heldProteinSets.ProteinDataHolder[i].namesSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]] << "\t - Pos : " << rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j] << "\t X: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]])) / 1000) << "\t Y: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]])) / 1000) << "\t Z: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 2 + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]])) / 1000) << std::endl;
									std::cout << "Atom B: " << heldProteinSets.ProteinDataHolder[i].namesSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]] << "\t - Pos : " << rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j] << "\t X: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]])) / 1000) << "\t Y: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]])) / 1000) << "\t Z: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 2 + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]])) / 1000) << std::endl << std::endl;

									std::cout << std::endl;
								}
								else if (outputType == 4)
								{
									filePrinter.printLineToOutputFile("Atom A: ", heldProteinSets.ProteinDataHolder[i].namesSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]], "\t - Pos : ", rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j], "\t X: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]])) / 1000), "\t Y: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]])) / 1000), "\t Z: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 2 + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]])) / 1000));
									filePrinter.printLineToOutputFile("Atom B: ", heldProteinSets.ProteinDataHolder[i].namesSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]], "\t - Pos : ", rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j], "\t X: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]])) / 1000), "\t Y: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]])) / 1000), "\t Z: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 2 + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]])) / 1000));
									filePrinter.printLineToOutputFile("");
								}

							}
						}

					}
					//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


				}

				//swap the pointer to which set needs to be processed and loaded next - so that the previously loaded set will be processed and a new set will be loaded into the spot of the already loaded set.
				switchLoadingAndProcessingSets(SetBeingProcessed, SetBeingLoaded);


			}
			m = clock();
			print_elapsed(n, m, "run time for bruteset mini: ");
			std::cout << std::endl;
		}

		filePrinter.closeOpenFile();
		
		
	}











	for (int i = 0; i < 5; i++) //For each of the 5 range lengths of stored protein: Unpin memory.
	{
		currentHeldEntries = heldProteinSets.ProteinDataHolder[i].heldEntries;
		if (currentHeldEntries > 0)
		{
			gpuErrchk(cudaHostUnregister(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets));
			gpuErrchk(cudaHostUnregister(heldProteinSets.ProteinDataHolder[i].kdTrees));
			gpuErrchk(cudaHostUnregister(heldProteinSets.ProteinDataHolder[i].namesSets));
		}
	}


	for (int i = 0; i < 2; i++)
	{
		gpuErrchk(cudaFree(rangeSearchSlots[i].d_atomACount));
		gpuErrchk(cudaFree(rangeSearchSlots[i].d_atomAPositionList));
		gpuErrchk(cudaFree(rangeSearchSlots[i].d_atomACurrentSearchDimensions));
		gpuErrchk(cudaFree(rangeSearchSlots[i].d_atomACurrentSearchKdTreePositions));
		gpuErrchk(cudaFree(rangeSearchSlots[i].d_atomAMatches));
		gpuErrchk(cudaFree(rangeSearchSlots[i].d_atomBMatches));
		gpuErrchk(cudaFree(rangeSearchSlots[i].d_MatchesCount));
		gpuErrchk(cudaFree(rangeSearchSlots[i].d_nextSearchCount));
		gpuErrchk(cudaFree(rangeSearchSlots[i].d_completionFlag));
		gpuErrchk(cudaFree(rangeSearchSlots[i].d_kdTreeSets));
		gpuErrchk(cudaFree(rangeSearchSlots[i].d_namesSets));
		gpuErrchk(cudaFree(rangeSearchSlots[i].d_xyzCoordsSets));
	}


	cudaStreamDestroy(streams[0]);
	cudaStreamDestroy(streams[1]);
	cudaDeviceReset();


};



void kdTreeSearchAllAvailableFilesSingleGpuLoadedEntryAtATime(rangeSearchSettings& settings, ProteinDataHandler heldProteinSets, AtomToNumHashTable atomReferenceTable)
{
	//all data has been preloaded into the host side ProteinDataHandler object. All that is needed is to initialise containers, moves sets of data to the gpu, process those sets of data and then return the results.
	std::cout << "PERFORMING TYPE 7 RANGE SEARCH: INDIVIDUAL LOAD GPU KD TREE" << std::endl;

	cudaStream_t streams[2];
	cudaStreamCreate(&streams[0]);
	cudaStreamCreate(&streams[1]);

	//Initialise host and device data holders
	int SetBeingProcessed = 0;
	int SetBeingLoaded = 1;
	gpuRangeSearchResources rangeSearchSlots[2];




	int IndividualEntryHolderSize = 16390 * 4; //The memory required to hold a max size coordinate array + 6 more atoms


	int blocks = 1;
	int threads = 1;
//	int currentEntry = 0;
	for (int i = 0; i < 2; i++)
	{
		gpuErrchk(cudaMalloc((void**)&rangeSearchSlots[i].d_atomACount, sizeof(int)));
		gpuErrchk(cudaMalloc((void**)&rangeSearchSlots[i].d_atomAPositionList, IndividualEntryHolderSize * 5));
		gpuErrchk(cudaMalloc((void**)&rangeSearchSlots[i].d_atomACurrentSearchDimensions, IndividualEntryHolderSize * 5));
		gpuErrchk(cudaMalloc((void**)&rangeSearchSlots[i].d_atomACurrentSearchKdTreePositions, IndividualEntryHolderSize * 5));
		gpuErrchk(cudaMalloc((void**)&rangeSearchSlots[i].d_atomAMatches, IndividualEntryHolderSize * 5));
		gpuErrchk(cudaMalloc((void**)&rangeSearchSlots[i].d_atomBMatches, IndividualEntryHolderSize * 5));
		gpuErrchk(cudaMalloc((void**)&rangeSearchSlots[i].d_MatchesCount, sizeof(int)));
		gpuErrchk(cudaMalloc((void**)&rangeSearchSlots[i].d_nextSearchCount, sizeof(int)));
		gpuErrchk(cudaMalloc((void**)&rangeSearchSlots[i].d_completionFlag, sizeof(int)));
		gpuErrchk(cudaMalloc((void**)&rangeSearchSlots[i].d_kdTreeSets, IndividualEntryHolderSize * 2));
		gpuErrchk(cudaMalloc((void**)&rangeSearchSlots[i].d_namesSets, IndividualEntryHolderSize / 2));
		gpuErrchk(cudaMalloc((void**)&rangeSearchSlots[i].d_xyzCoordsSets, IndividualEntryHolderSize * 3));
		rangeSearchSlots[i].h_atomACount = (int*)malloc(sizeof(int));
		rangeSearchSlots[i].h_atomAPositionList = (int*)malloc(IndividualEntryHolderSize * 5);
		rangeSearchSlots[i].h_nextSearchCount = (int*)malloc(sizeof(int));
		rangeSearchSlots[i].h_completionFlag = (int*)malloc(sizeof(int));
		rangeSearchSlots[i].h_MatchesCount = (int*)malloc(sizeof(int));
		rangeSearchSlots[i].h_atomAMatches = (int*)malloc(IndividualEntryHolderSize * 5);
		rangeSearchSlots[i].h_atomBMatches = (int*)malloc(IndividualEntryHolderSize * 5);
	}




	int soughtAtomANumber = atomReferenceTable.retrieveHashValue(settings.AtomTypeOne);
	int soughtAtomBNumber = atomReferenceTable.retrieveHashValue(settings.AtomTypeTwo);

	int currentMaxEntrySize;
	int currentHeldEntries;


	for (int i = 0; i < 5; i++) //For each of the 5 range lengths of stored protein: Pin memory.
	{
		currentMaxEntrySize = heldProteinSets.ProteinDataHolder[i].MaxEntrySize;
		currentHeldEntries = heldProteinSets.ProteinDataHolder[i].heldEntries;
		if (currentHeldEntries > 0)
		{
			gpuErrchk(cudaHostRegister(heldProteinSets.ProteinDataHolder[i].kdTrees, currentMaxEntrySize * currentHeldEntries * sizeof(int) * 2, 0));
			gpuErrchk(cudaHostRegister(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets, currentMaxEntrySize * currentHeldEntries * sizeof(int) * 3, 0));
			gpuErrchk(cudaHostRegister(heldProteinSets.ProteinDataHolder[i].namesSets, currentMaxEntrySize * currentHeldEntries * sizeof(short), 0));
		}
	}



	for (int i = 0; i < 5; i++)//For each of the 5 range lengths of stored protein:
	{
		currentMaxEntrySize = heldProteinSets.ProteinDataHolder[i].MaxEntrySize;
		int TotalEntries = heldProteinSets.ProteinDataHolder[i].heldEntries;

		clock_t n, m;



		int outputType = settings.resultsPrintFormat;
		outputHandler filePrinter;
		std::string printType;
		if (outputType == 3) { printType = "_Summary"; }
		else if (outputType == 4) { printType = "_Detailed"; }

		if (outputType == 3 || outputType == 4)	{ filePrinter.initializeOutputfile("GpuKdTreeResults_Range_", currentMaxEntrySize, "_Files_", TotalEntries, printType); }

		n = clock();

		if (TotalEntries > 0)
		{
			std::cout << "Processing Range set: " << i << std::endl;
			std::cout << "Number of present entries is: " << heldProteinSets.ProteinDataHolder[i].heldEntries << std::endl;
			for (int currentEntry = 0; currentEntry < TotalEntries + 1; currentEntry++)
			{
				if (currentEntry == 0)
				{

					//Load first set of details onto the gpu but do not process them 

					//gpuErrchk(cudaMemcpy(rangeSearchSlots[SetBeingLoaded].d_xyzCoordsSets, heldProteinSets.ProteinDataHolder[i].xyzCoordsSets + currentMaxEntrySize*currentEntry * 3, currentMaxEntrySize*sizeof(int) * 3, cudaMemcpyHostToDevice));
					gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingLoaded].d_xyzCoordsSets, heldProteinSets.ProteinDataHolder[i].xyzCoordsSets + currentMaxEntrySize*currentEntry * 3, currentMaxEntrySize*sizeof(int) * 3, cudaMemcpyHostToDevice, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingLoaded].d_namesSets, heldProteinSets.ProteinDataHolder[i].namesSets + currentMaxEntrySize*currentEntry, currentMaxEntrySize *sizeof(short), cudaMemcpyHostToDevice, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingLoaded].d_kdTreeSets, heldProteinSets.ProteinDataHolder[i].kdTrees + currentMaxEntrySize*currentEntry * 2, currentMaxEntrySize * 2 * sizeof(int), cudaMemcpyHostToDevice, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_atomACount, 0, sizeof(int), streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_atomAPositionList, -1, IndividualEntryHolderSize * 5, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_atomAMatches, -1, IndividualEntryHolderSize * 5, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_atomBMatches, -1, IndividualEntryHolderSize * 5, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_atomACurrentSearchDimensions, 0, IndividualEntryHolderSize * 5, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_atomACurrentSearchKdTreePositions, -1, IndividualEntryHolderSize * 5, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_MatchesCount, 0, sizeof(int), streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_nextSearchCount, 0, sizeof(int), streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_completionFlag, 0, sizeof(int), streams[SetBeingLoaded]));
					rangeSearchSlots[SetBeingLoaded].h_MatchesCount[0] = 0;




				}
				else if (currentEntry == TotalEntries)
				{
					//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					//process current loaded set
					threads = 512;
					blocks = heldProteinSets.ProteinDataHolder[i].MaxEntrySize / 512;
					device_side__locateElement << <blocks, threads, 0, streams[SetBeingProcessed] >> >(rangeSearchSlots[SetBeingProcessed].d_namesSets, 0, soughtAtomANumber, rangeSearchSlots[SetBeingProcessed].d_atomAPositionList, rangeSearchSlots[SetBeingProcessed].d_atomACount, rangeSearchSlots[SetBeingProcessed].d_atomACurrentSearchDimensions, rangeSearchSlots[SetBeingProcessed].d_atomACurrentSearchKdTreePositions, heldProteinSets.ProteinDataHolder[i].MaxEntrySize);

					gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingProcessed].h_atomACount, rangeSearchSlots[SetBeingProcessed].d_atomACount, sizeof(int), cudaMemcpyDeviceToHost, streams[SetBeingProcessed]));
					//gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingProcessed].h_atomAPositionList, rangeSearchSlots[SetBeingProcessed].d_atomAPositionList, safeHolderSize, cudaMemcpyDeviceToHost, streams[SetBeingProcessed]));

					int maxNumberoFCycles = 13 + (heldProteinSets.ProteinDataHolder[i].MaxEntrySize) / 1024;
					rangeSearchSlots[SetBeingProcessed].blocks = 1;
					rangeSearchSlots[SetBeingProcessed].threads = 1;
					if (rangeSearchSlots[SetBeingProcessed].h_atomACount[0] > 0)
					{
						for (int p = 0; p < maxNumberoFCycles; p++)
						{

							gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingProcessed].h_atomACount, rangeSearchSlots[SetBeingProcessed].d_atomACount, sizeof(int), cudaMemcpyDeviceToHost, streams[SetBeingProcessed]));

							if (rangeSearchSlots[SetBeingProcessed].h_atomACount[0] < 513)
							{
								rangeSearchSlots[SetBeingProcessed].blocks = 1;
								rangeSearchSlots[SetBeingProcessed].threads = rangeSearchSlots[SetBeingProcessed].h_atomACount[0];
							}
							else
							{
								rangeSearchSlots[SetBeingProcessed].threads = 512;
								rangeSearchSlots[SetBeingProcessed].blocks = rangeSearchSlots[SetBeingProcessed].h_atomACount[0] / 512 + 1;
							}



							///seems to have redundant entries
							//gpuErrchk(cudaMemcpy(rangeSearchSlots[SetBeingProcessed].h_MatchesCount, rangeSearchSlots[SetBeingProcessed].d_MatchesCount, sizeof(int), cudaMemcpyDeviceToHost));//For testing...
							device_side_ProcessCurrentTreePositionsV3 << <rangeSearchSlots[SetBeingProcessed].blocks, rangeSearchSlots[SetBeingProcessed].threads, 0, streams[SetBeingProcessed] >> >(rangeSearchSlots[SetBeingProcessed].d_xyzCoordsSets, rangeSearchSlots[SetBeingProcessed].d_namesSets, rangeSearchSlots[SetBeingProcessed].d_kdTreeSets, rangeSearchSlots[SetBeingProcessed].d_atomAMatches, rangeSearchSlots[SetBeingProcessed].d_atomBMatches, rangeSearchSlots[SetBeingProcessed].d_MatchesCount, rangeSearchSlots[SetBeingProcessed].d_atomACurrentSearchKdTreePositions, rangeSearchSlots[SetBeingProcessed].d_atomAPositionList, rangeSearchSlots[SetBeingProcessed].d_atomACurrentSearchDimensions, rangeSearchSlots[SetBeingProcessed].d_atomACount, heldProteinSets.ProteinDataHolder[i].MaxEntrySize, settings.requiredProximity*settings.requiredProximity, soughtAtomBNumber, rangeSearchSlots[SetBeingProcessed].d_nextSearchCount, 0, heldProteinSets.ProteinDataHolder[i].KdTreeSize, settings.requiredProximity, p);





							//gpuErrchk(cudaMemcpy(rangeSearch.h_MatchesCount, rangeSearch.d_MatchesCount, sizeof(int), cudaMemcpyDeviceToHost));//For testing...
							SetNextCountAsCurrentAndCheckFlag << <1, 1, 0, streams[SetBeingProcessed] >> >(rangeSearchSlots[SetBeingProcessed].d_atomACount, rangeSearchSlots[SetBeingProcessed].d_nextSearchCount, rangeSearchSlots[SetBeingProcessed].d_completionFlag); //update device side counters reflecting how many search items exist in the list


							gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingProcessed].h_completionFlag, rangeSearchSlots[SetBeingProcessed].d_completionFlag, sizeof(int), cudaMemcpyDeviceToHost, streams[SetBeingProcessed]));
							if (rangeSearchSlots[SetBeingProcessed].h_completionFlag[0] == 1)
							{
								p = 1000;
							}

						}
					}

					gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingProcessed].h_MatchesCount, rangeSearchSlots[SetBeingProcessed].d_MatchesCount, sizeof(int), cudaMemcpyDeviceToHost, streams[SetBeingProcessed]));

					if (outputType == 1 || outputType == 2)
						std::cout << "Number of matches in file " << currentEntry << " in set " << i << " is: " << rangeSearchSlots[SetBeingProcessed].h_MatchesCount[0] << std::endl;
					else if (outputType == 3 || outputType == 4)
						filePrinter.printLineToOutputFile("Number of matches in file ", currentEntry, " in set ", i, "  is: ", rangeSearchSlots[SetBeingProcessed].h_MatchesCount[0]);







					if ((rangeSearchSlots[SetBeingProcessed].h_MatchesCount[0] > 0) && (rangeSearchSlots[SetBeingProcessed].h_atomACount[0] > 0))
					{
						gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingProcessed].h_atomAMatches, rangeSearchSlots[SetBeingProcessed].d_atomAMatches, sizeof(int)*rangeSearchSlots[SetBeingProcessed].h_MatchesCount[0], cudaMemcpyDeviceToHost, streams[SetBeingProcessed]));
						gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingProcessed].h_atomBMatches, rangeSearchSlots[SetBeingProcessed].d_atomBMatches, sizeof(int)*rangeSearchSlots[SetBeingProcessed].h_MatchesCount[0], cudaMemcpyDeviceToHost, streams[SetBeingProcessed]));
						int valuesOffset = heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 3 * currentEntry;

						for (int j = 0; j < rangeSearchSlots[SetBeingProcessed].h_MatchesCount[0]; j++)
						{
							if (outputType == 2)
							{
								std::cout << "Atom A: " << heldProteinSets.ProteinDataHolder[i].namesSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]] << "\t - Pos : " << rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j] << "\t X: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]])) / 1000) << "\t Y: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]])) / 1000) << "\t Z: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 2 + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]])) / 1000) << std::endl;
								std::cout << "Atom B: " << heldProteinSets.ProteinDataHolder[i].namesSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]] << "\t - Pos : " << rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j] << "\t X: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]])) / 1000) << "\t Y: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]])) / 1000) << "\t Z: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 2 + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]])) / 1000) << std::endl << std::endl;

								std::cout << std::endl;
							}
							else if (outputType == 4)
							{
								filePrinter.printLineToOutputFile("Atom A: ", heldProteinSets.ProteinDataHolder[i].namesSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]], "\t - Pos : ", rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j], "\t X: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]])) / 1000), "\t Y: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]])) / 1000), "\t Z: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 2 + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]])) / 1000));
								filePrinter.printLineToOutputFile("Atom B: ", heldProteinSets.ProteinDataHolder[i].namesSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]], "\t - Pos : ", rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j], "\t X: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]])) / 1000), "\t Y: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]])) / 1000), "\t Z: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 2 + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]])) / 1000));
								filePrinter.printLineToOutputFile("");
							}

						}
					}
					//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


				}
				else
				{


					//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					//Load next set

					gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingLoaded].d_namesSets, heldProteinSets.ProteinDataHolder[i].namesSets + currentMaxEntrySize*currentEntry, currentMaxEntrySize *sizeof(short), cudaMemcpyHostToDevice, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingLoaded].d_kdTreeSets, heldProteinSets.ProteinDataHolder[i].kdTrees + currentMaxEntrySize*currentEntry * 2, currentMaxEntrySize * 2 * sizeof(int), cudaMemcpyHostToDevice, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_atomACount, 0, sizeof(int), streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_atomAPositionList, -1, IndividualEntryHolderSize * 5, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_atomAMatches, -1, IndividualEntryHolderSize * 5, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_atomBMatches, -1, IndividualEntryHolderSize * 5, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_atomACurrentSearchDimensions, 0, IndividualEntryHolderSize * 5, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_atomACurrentSearchKdTreePositions, -1, IndividualEntryHolderSize * 5, streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_MatchesCount, 0, sizeof(int), streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_nextSearchCount, 0, sizeof(int), streams[SetBeingLoaded]));
					gpuErrchk(cudaMemsetAsync(rangeSearchSlots[SetBeingLoaded].d_completionFlag, 0, sizeof(int), streams[SetBeingLoaded]));
					rangeSearchSlots[SetBeingProcessed].h_MatchesCount[0] = 0;
					gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingLoaded].d_xyzCoordsSets, heldProteinSets.ProteinDataHolder[i].xyzCoordsSets + currentMaxEntrySize*currentEntry * 3, currentMaxEntrySize*sizeof(int) * 3, cudaMemcpyHostToDevice, streams[SetBeingLoaded]));


					//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					//process current loaded set
					threads = 512;
					blocks = heldProteinSets.ProteinDataHolder[i].MaxEntrySize / 512;
					device_side__locateElement << <blocks, threads, 0, streams[SetBeingProcessed] >> >(rangeSearchSlots[SetBeingProcessed].d_namesSets, 0, soughtAtomANumber, rangeSearchSlots[SetBeingProcessed].d_atomAPositionList, rangeSearchSlots[SetBeingProcessed].d_atomACount, rangeSearchSlots[SetBeingProcessed].d_atomACurrentSearchDimensions, rangeSearchSlots[SetBeingProcessed].d_atomACurrentSearchKdTreePositions, heldProteinSets.ProteinDataHolder[i].MaxEntrySize);

					gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingProcessed].h_atomACount, rangeSearchSlots[SetBeingProcessed].d_atomACount, sizeof(int), cudaMemcpyDeviceToHost, streams[SetBeingProcessed]));
					//gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingProcessed].h_atomAPositionList, rangeSearchSlots[SetBeingProcessed].d_atomAPositionList, safeHolderSize, cudaMemcpyDeviceToHost, streams[SetBeingProcessed]));

					int maxNumberoFCycles = 13 + (heldProteinSets.ProteinDataHolder[i].MaxEntrySize) / 1024;
					rangeSearchSlots[SetBeingProcessed].blocks = 1;
					rangeSearchSlots[SetBeingProcessed].threads = 1;
					if (rangeSearchSlots[SetBeingProcessed].h_atomACount[0] > 0)
					{
						for (int p = 0; p < maxNumberoFCycles; p++)
						{

							gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingProcessed].h_atomACount, rangeSearchSlots[SetBeingProcessed].d_atomACount, sizeof(int), cudaMemcpyDeviceToHost, streams[SetBeingProcessed]));

							if (rangeSearchSlots[SetBeingProcessed].h_atomACount[0] < 513)
							{
								rangeSearchSlots[SetBeingProcessed].blocks = 1;
								rangeSearchSlots[SetBeingProcessed].threads = rangeSearchSlots[SetBeingProcessed].h_atomACount[0];
							}
							else
							{
								rangeSearchSlots[SetBeingProcessed].threads = 512;
								rangeSearchSlots[SetBeingProcessed].blocks = rangeSearchSlots[SetBeingProcessed].h_atomACount[0] / 512 + 1;
							}




							//gpuErrchk(cudaMemcpy(rangeSearchSlots[SetBeingProcessed].h_MatchesCount, rangeSearchSlots[SetBeingProcessed].d_MatchesCount, sizeof(int), cudaMemcpyDeviceToHost));//For testing...
							device_side_ProcessCurrentTreePositionsV3 << <rangeSearchSlots[SetBeingProcessed].blocks, rangeSearchSlots[SetBeingProcessed].threads, 0, streams[SetBeingProcessed] >> >(rangeSearchSlots[SetBeingProcessed].d_xyzCoordsSets, rangeSearchSlots[SetBeingProcessed].d_namesSets, rangeSearchSlots[SetBeingProcessed].d_kdTreeSets, rangeSearchSlots[SetBeingProcessed].d_atomAMatches, rangeSearchSlots[SetBeingProcessed].d_atomBMatches, rangeSearchSlots[SetBeingProcessed].d_MatchesCount, rangeSearchSlots[SetBeingProcessed].d_atomACurrentSearchKdTreePositions, rangeSearchSlots[SetBeingProcessed].d_atomAPositionList, rangeSearchSlots[SetBeingProcessed].d_atomACurrentSearchDimensions, rangeSearchSlots[SetBeingProcessed].d_atomACount, heldProteinSets.ProteinDataHolder[i].MaxEntrySize, settings.requiredProximity*settings.requiredProximity, soughtAtomBNumber, rangeSearchSlots[SetBeingProcessed].d_nextSearchCount, 0, heldProteinSets.ProteinDataHolder[i].KdTreeSize, settings.requiredProximity, p);





							//gpuErrchk(cudaMemcpy(rangeSearch.h_MatchesCount, rangeSearch.d_MatchesCount, sizeof(int), cudaMemcpyDeviceToHost));//For testing...
							SetNextCountAsCurrentAndCheckFlag << <1, 1, 0, streams[SetBeingProcessed] >> >(rangeSearchSlots[SetBeingProcessed].d_atomACount, rangeSearchSlots[SetBeingProcessed].d_nextSearchCount, rangeSearchSlots[SetBeingProcessed].d_completionFlag); //update device side counters reflecting how many search items exist in the list


							gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingProcessed].h_completionFlag, rangeSearchSlots[SetBeingProcessed].d_completionFlag, sizeof(int), cudaMemcpyDeviceToHost, streams[SetBeingProcessed]));
							if (rangeSearchSlots[SetBeingProcessed].h_completionFlag[0] == 1)
							{
								p = 1000;
							}

						}
					}

					gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingProcessed].h_MatchesCount, rangeSearchSlots[SetBeingProcessed].d_MatchesCount, sizeof(int), cudaMemcpyDeviceToHost, streams[SetBeingProcessed]));

					if (outputType == 1 || outputType == 2)
						std::cout << "Number of matches in file " << currentEntry << " in set " << i << " is: " << rangeSearchSlots[SetBeingProcessed].h_MatchesCount[0] << std::endl;
					else if (outputType == 3 || outputType == 4)
						filePrinter.printLineToOutputFile("Number of matches in file ", currentEntry, " in set ", i, "  is: ", rangeSearchSlots[SetBeingProcessed].h_MatchesCount[0]);







					if ((rangeSearchSlots[SetBeingProcessed].h_MatchesCount[0] > 0) && (rangeSearchSlots[SetBeingProcessed].h_atomACount[0] > 0))
					{
						gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingProcessed].h_atomAMatches, rangeSearchSlots[SetBeingProcessed].d_atomAMatches, sizeof(int)*rangeSearchSlots[SetBeingProcessed].h_MatchesCount[0], cudaMemcpyDeviceToHost, streams[SetBeingProcessed]));
						gpuErrchk(cudaMemcpyAsync(rangeSearchSlots[SetBeingProcessed].h_atomBMatches, rangeSearchSlots[SetBeingProcessed].d_atomBMatches, sizeof(int)*rangeSearchSlots[SetBeingProcessed].h_MatchesCount[0], cudaMemcpyDeviceToHost, streams[SetBeingProcessed]));
						int valuesOffset = heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 3 * currentEntry;

						for (int j = 0; j < rangeSearchSlots[SetBeingProcessed].h_MatchesCount[0]; j++)
						{
							if (outputType == 2)
							{
								std::cout << "Atom A: " << heldProteinSets.ProteinDataHolder[i].namesSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]] << "\t - Pos : " << rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j] << "\t X: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]])) / 1000) << "\t Y: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]])) / 1000) << "\t Z: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 2 + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]])) / 1000) << std::endl;
								std::cout << "Atom B: " << heldProteinSets.ProteinDataHolder[i].namesSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]] << "\t - Pos : " << rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j] << "\t X: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]])) / 1000) << "\t Y: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]])) / 1000) << "\t Z: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 2 + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]])) / 1000) << std::endl << std::endl;

								std::cout << std::endl;
							}
							else if (outputType == 4)
							{
								filePrinter.printLineToOutputFile("Atom A: ", heldProteinSets.ProteinDataHolder[i].namesSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]], "\t - Pos : ", rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j], "\t X: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]])) / 1000), "\t Y: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]])) / 1000), "\t Z: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 2 + rangeSearchSlots[SetBeingProcessed].h_atomAMatches[j]])) / 1000));
								filePrinter.printLineToOutputFile("Atom B: ", heldProteinSets.ProteinDataHolder[i].namesSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]], "\t - Pos : ", rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j], "\t X: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]])) / 1000), "\t Y: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]])) / 1000), "\t Z: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 2 + rangeSearchSlots[SetBeingProcessed].h_atomBMatches[j]])) / 1000));
								filePrinter.printLineToOutputFile("");
							}

						}
					}
					//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


				}

				//swap the pointer to which set needs to be processed and loaded next - so that the previously loaded set will be processed and a new set will be loaded into the spot of the already loaded set.
				switchLoadingAndProcessingSets(SetBeingProcessed, SetBeingLoaded);



			}
			m = clock();
			print_elapsed(n, m, "run time for bruteset mini: ");
			std::cout << std::endl;
		}

		filePrinter.closeOpenFile();
		

	}




	for (int i = 0; i < 5; i++) //For each of the 5 range lengths of stored protein: Unpin memory.
	{
		currentHeldEntries = heldProteinSets.ProteinDataHolder[i].heldEntries;
		if (currentHeldEntries > 0)
		{
			gpuErrchk(cudaHostUnregister(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets));
			gpuErrchk(cudaHostUnregister(heldProteinSets.ProteinDataHolder[i].kdTrees));
			gpuErrchk(cudaHostUnregister(heldProteinSets.ProteinDataHolder[i].namesSets));
		}
	}


	for (int i = 0; i < 2; i++)
	{
		gpuErrchk(cudaFree(rangeSearchSlots[i].d_atomACount));
		gpuErrchk(cudaFree(rangeSearchSlots[i].d_atomAPositionList));
		gpuErrchk(cudaFree(rangeSearchSlots[i].d_atomACurrentSearchDimensions));
		gpuErrchk(cudaFree(rangeSearchSlots[i].d_atomACurrentSearchKdTreePositions));
		gpuErrchk(cudaFree(rangeSearchSlots[i].d_atomAMatches));
		gpuErrchk(cudaFree(rangeSearchSlots[i].d_atomBMatches));
		gpuErrchk(cudaFree(rangeSearchSlots[i].d_MatchesCount));
		gpuErrchk(cudaFree(rangeSearchSlots[i].d_nextSearchCount));
		gpuErrchk(cudaFree(rangeSearchSlots[i].d_completionFlag));
		gpuErrchk(cudaFree(rangeSearchSlots[i].d_kdTreeSets));
		gpuErrchk(cudaFree(rangeSearchSlots[i].d_namesSets));
		gpuErrchk(cudaFree(rangeSearchSlots[i].d_xyzCoordsSets));
	}


	cudaStreamDestroy(streams[0]);
	cudaStreamDestroy(streams[1]);


}

//for small sets
void gpuKdTreeUnoptimisedRangeSearchAllLoadedSets(rangeSearchSettings &settings, ProteinDataHandler heldProteinSets, AtomToNumHashTable atomReferenceTable)
{
	//all data has been preloaded into the host side ProteinDataHandler object. All that is needed is to initialise containers, moves sets of data to the gpu, process those sets of data and then return the results.
	std::cout << "PERFORMING TYPE 6 RANGE SEARCH: SINGLE BULK BATCH GPU KD-TREE" << std::endl;

	//Initialise host and device data holders

	gpuRangeSearchResources rangeSearch;
	int safeHolderSize = 32000 * 90 * sizeof(int) * 3; //This will crash eventually - cant fit massive sets on the gpu in one go -.-

	int blocks = 1;
	int threads = 1;
//	int currentEntry = 0;
	gpuErrchk(cudaMalloc((void**)&rangeSearch.d_atomACount, sizeof(int)));
	gpuErrchk(cudaMalloc((void**)&rangeSearch.d_atomAPositionList, safeHolderSize));
	gpuErrchk(cudaMalloc((void**)&rangeSearch.d_atomACurrentSearchDimensions, safeHolderSize));
	gpuErrchk(cudaMalloc((void**)&rangeSearch.d_atomACurrentSearchKdTreePositions, safeHolderSize));
	gpuErrchk(cudaMalloc((void**)&rangeSearch.d_atomAMatches, safeHolderSize));
	gpuErrchk(cudaMalloc((void**)&rangeSearch.d_atomBMatches, safeHolderSize));
	gpuErrchk(cudaMalloc((void**)&rangeSearch.d_MatchesCount, sizeof(int)));
	gpuErrchk(cudaMalloc((void**)&rangeSearch.d_nextSearchCount, sizeof(int)));
	gpuErrchk(cudaMalloc((void**)&rangeSearch.d_completionFlag, sizeof(int)));
	gpuErrchk(cudaMalloc((void**)&rangeSearch.d_kdTreeSets, safeHolderSize * 2));
	gpuErrchk(cudaMalloc((void**)&rangeSearch.d_namesSets, safeHolderSize / sizeof(int) * sizeof(short)));
	gpuErrchk(cudaMalloc((void**)&rangeSearch.d_xyzCoordsSets, safeHolderSize * 3));
	rangeSearch.h_atomACount = (int*)malloc(sizeof(int));
	//rangeSearch.h_atomACurrentSearchDimensions = (int*)malloc(safeHolderSize);
	//rangeSearch.h_atomACurrentSearchKdTreePositions = (int*)malloc(safeHolderSize);
	rangeSearch.h_atomAPositionList = (int*)malloc(safeHolderSize);
	rangeSearch.h_nextSearchCount = (int*)malloc(sizeof(int));
	rangeSearch.h_completionFlag = (int*)malloc(sizeof(int));
	rangeSearch.h_MatchesCount = (int*)malloc(sizeof(int));
	rangeSearch.h_atomAMatches = (int*)malloc(safeHolderSize);
	rangeSearch.h_atomBMatches = (int*)malloc(safeHolderSize);

	cudaStream_t streams[3];

	cudaStreamCreate(&streams[0]);
	cudaStreamCreate(&streams[1]);
	cudaStreamCreate(&streams[2]);

	//load data into holder arrays:
	int entrySize;
	int kdEntrySize;
	int proteinsInSet;


	//initiate search


	//int soughtAtomANumber = 1414;
	//int soughtAtomBNumber = 1514;

	int soughtAtomANumber = atomReferenceTable.retrieveHashValue(settings.AtomTypeOne);//heldProteinSets.ProteinDataHolder[0].namesSets[0];//
	int soughtAtomBNumber = atomReferenceTable.retrieveHashValue(settings.AtomTypeTwo);//heldProteinSets.ProteinDataHolder[0].namesSets[1];//
	int maxDistanceSquared = settings.requiredProximity*settings.requiredProximity;


	for (int i = 0; i < 5; i++)
	{
		entrySize = heldProteinSets.ProteinDataHolder[i].MaxEntrySize;
		kdEntrySize = heldProteinSets.ProteinDataHolder[i].KdTreeSize;
		proteinsInSet = heldProteinSets.ProteinDataHolder[i].heldEntries;
		
		if (proteinsInSet > 0)
		{
			gpuErrchk(cudaHostRegister(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets, (proteinsInSet*entrySize*sizeof(int) * 3), 0));
			gpuErrchk(cudaHostRegister(heldProteinSets.ProteinDataHolder[i].namesSets, proteinsInSet*entrySize *sizeof(short), 0));
			gpuErrchk(cudaHostRegister(heldProteinSets.ProteinDataHolder[i].kdTrees, proteinsInSet*kdEntrySize*sizeof(int), 0));
		}
	}



	for (int i = 0; i < 5; i++)
	{
		entrySize = heldProteinSets.ProteinDataHolder[i].MaxEntrySize;
		kdEntrySize = heldProteinSets.ProteinDataHolder[i].KdTreeSize;
		proteinsInSet = heldProteinSets.ProteinDataHolder[i].heldEntries;
		//the x,y and z data is loaded from the 3 host arrays into a single gpu array - the xyz array
		//gpuErrchk(cudaMemcpyAsync(rangeSearch.d_xyzCoordsSets, heldProteinSets.ProteinDataHolder[i].xCoordsSets, entrySize*sizeof(int), cudaMemcpyHostToDevice/*, streams[0]*/)); //seeing as the entire x lot is put before the entire y lot, i could probably still use this approach.
		//gpuErrchk(cudaMemcpyAsync(rangeSearch.d_xyzCoordsSets + entrySize, heldProteinSets.ProteinDataHolder[i].yCoordsSets, entrySize*sizeof(int), cudaMemcpyHostToDevice/*, streams[1]*/));
		//gpuErrchk(cudaMemcpyAsync(rangeSearch.d_xyzCoordsSets + entrySize * 2, heldProteinSets.ProteinDataHolder[i].zCoordsSets, entrySize*sizeof(int), cudaMemcpyHostToDevice/*, streams[2]*/));
		if (proteinsInSet > 0)
		{
			gpuErrchk(cudaMemcpy(rangeSearch.d_xyzCoordsSets, heldProteinSets.ProteinDataHolder[i].xyzCoordsSets, proteinsInSet*entrySize*sizeof(int) * 3, cudaMemcpyHostToDevice));
			gpuErrchk(cudaMemcpyAsync(rangeSearch.d_namesSets, heldProteinSets.ProteinDataHolder[i].namesSets, proteinsInSet*entrySize *sizeof(short), cudaMemcpyHostToDevice/*, streams[0]*/));
			gpuErrchk(cudaMemcpyAsync(rangeSearch.d_kdTreeSets, heldProteinSets.ProteinDataHolder[i].kdTrees, proteinsInSet*kdEntrySize*sizeof(int), cudaMemcpyHostToDevice/*, streams[0]*/));
			std::cout << "Processing Range set: " << i << std::endl;
			std::cout << "Number of present entries is: " << proteinsInSet << std::endl;


			int outputType = settings.resultsPrintFormat;
			outputHandler filePrinter;
			std::string printType;

			if (outputType == 3)
			{
				printType = "_Summary";
			}
			else if (outputType == 4)
			{
				printType = "_Detailed";
			}


			if (outputType == 3 || outputType == 4)//will move into loop shortly
			{
				filePrinter.initializeOutputfile("GpuKdResults_Range_", entrySize, "_Files_", proteinsInSet, printType);
			}







			for (int currentEntry = 0; currentEntry < heldProteinSets.ProteinDataHolder[i].heldEntries; currentEntry++)
			{
				gpuErrchk(cudaMemsetAsync(rangeSearch.d_atomACount, 0, sizeof(int)/*, streams[0]*/));
				gpuErrchk(cudaMemsetAsync(rangeSearch.d_atomAPositionList, -1, safeHolderSize/*, streams[1]*/));
				gpuErrchk(cudaMemsetAsync(rangeSearch.d_atomAMatches, -1, safeHolderSize/*, streams[2]*/));
				gpuErrchk(cudaMemsetAsync(rangeSearch.d_atomBMatches, -1, safeHolderSize/*, streams[0]*/));
				gpuErrchk(cudaMemsetAsync(rangeSearch.d_atomACurrentSearchDimensions, 0, safeHolderSize/*, streams[1]*/));
				gpuErrchk(cudaMemsetAsync(rangeSearch.d_atomACurrentSearchKdTreePositions, -1, safeHolderSize/*, streams[2]*/));
				gpuErrchk(cudaMemsetAsync(rangeSearch.d_MatchesCount, 0, sizeof(int)/*, streams[0]*/));
				gpuErrchk(cudaMemsetAsync(rangeSearch.d_nextSearchCount, 0, sizeof(int)/*, streams[1]*/));
				gpuErrchk(cudaMemsetAsync(rangeSearch.d_completionFlag, 0, sizeof(int)/*, streams[2]*/));
				rangeSearch.h_MatchesCount[0] = 0;


				cudaDeviceSynchronize();
				threads = 512;
				blocks = heldProteinSets.ProteinDataHolder[i].MaxEntrySize / 512;
				device_side__locateElement << <blocks, threads/*, 0, streams[0] */ >> >(rangeSearch.d_namesSets, currentEntry, soughtAtomANumber, rangeSearch.d_atomAPositionList, rangeSearch.d_atomACount, rangeSearch.d_atomACurrentSearchDimensions, rangeSearch.d_atomACurrentSearchKdTreePositions, heldProteinSets.ProteinDataHolder[i].MaxEntrySize);

				gpuErrchk(cudaMemcpyAsync(rangeSearch.h_atomACount, rangeSearch.d_atomACount, sizeof(int), cudaMemcpyDeviceToHost, streams[0]));
				gpuErrchk(cudaMemcpyAsync(rangeSearch.h_atomAPositionList, rangeSearch.d_atomAPositionList, safeHolderSize, cudaMemcpyDeviceToHost, streams[1]));
				// ---- std::cout << "File: " << currentEntry << " - atom A count: " << rangeSearch.h_atomACount[0] << "  ";

				int maxNumberoFCycles = 13 + (heldProteinSets.ProteinDataHolder[i].MaxEntrySize) / 1024;
				rangeSearch.blocks = 1;
				rangeSearch.threads = 1;
				if (rangeSearch.h_atomACount[0] > 0)
				{
					//std::cout << std::endl << "Entry " << currentEntry << " Contains : " << rangeSearch.h_atomACount[0] << " unique instances of element : " << soughtAtomANumber << std::endl;// If any searches are being performed on the kd tree, this prints out how many instances of the first atom are present in the tree.
					for (int p = 0; p < maxNumberoFCycles; p++)
					{
						//gpuErrchk(cudaMemcpy(rangeSearch.h_atomAPositionList, rangeSearch.d_atomAPositionList, sizeof(int) * rangeSearch.h_atomACount[0], cudaMemcpyDeviceToHost));// for testing only



						gpuErrchk(cudaMemcpy(rangeSearch.h_atomACount, rangeSearch.d_atomACount, sizeof(int) * 1, cudaMemcpyDeviceToHost));
						//std::cout << "Current active search count: " << rangeSearch.h_atomACount[0] << std::endl; //This prints how many searches are currently being proccessed on the current kd-tree

						//calculateInitialBlocksAndThreads(blocks, threads, heldProteinSets.ProteinDataHolder[i].MaxEntrySize);

						if (rangeSearch.h_atomACount[0] < 513)
						{
							rangeSearch.blocks = 1;
							rangeSearch.threads = rangeSearch.h_atomACount[0];
						}
						else
						{
							rangeSearch.threads = 512;
							rangeSearch.blocks = rangeSearch.h_atomACount[0] / 512 + 1;
						}


						gpuErrchk(cudaMemcpy(rangeSearch.h_MatchesCount, rangeSearch.d_MatchesCount, sizeof(int), cudaMemcpyDeviceToHost));//For testing...
						device_side_ProcessCurrentTreePositionsV3 << <rangeSearch.blocks, rangeSearch.threads >> >(rangeSearch.d_xyzCoordsSets, rangeSearch.d_namesSets, rangeSearch.d_kdTreeSets, rangeSearch.d_atomAMatches, rangeSearch.d_atomBMatches, rangeSearch.d_MatchesCount, rangeSearch.d_atomACurrentSearchKdTreePositions, rangeSearch.d_atomAPositionList, rangeSearch.d_atomACurrentSearchDimensions, rangeSearch.d_atomACount, heldProteinSets.ProteinDataHolder[i].MaxEntrySize, maxDistanceSquared, soughtAtomBNumber, rangeSearch.d_nextSearchCount, currentEntry, heldProteinSets.ProteinDataHolder[i].KdTreeSize, settings.requiredProximity, p);
						gpuErrchk(cudaMemcpy(rangeSearch.h_MatchesCount, rangeSearch.d_MatchesCount, sizeof(int), cudaMemcpyDeviceToHost));//For testing...
						SetNextCountAsCurrentAndCheckFlag << <1, 1 >> >(rangeSearch.d_atomACount, rangeSearch.d_nextSearchCount, rangeSearch.d_completionFlag); //update device side counters reflecting how many search items exist in the list


						gpuErrchk(cudaMemcpy(rangeSearch.h_completionFlag, rangeSearch.d_completionFlag, sizeof(int), cudaMemcpyDeviceToHost));
						if (rangeSearch.h_completionFlag[0] == 1)
						{
							p = 1000;
						}

					}
				}

				gpuErrchk(cudaMemcpy(rangeSearch.h_MatchesCount, rangeSearch.d_MatchesCount, sizeof(int), cudaMemcpyDeviceToHost));

				if (outputType == 1 || outputType == 2)
					std::cout << "Number of matches in file " << currentEntry << " in set " << i << " is: " << rangeSearch.h_MatchesCount[0] << std::endl;
				else if (outputType == 3 || outputType == 4)
					filePrinter.printLineToOutputFile("Number of matches in file ", currentEntry, " in set ", i, "  is: ", rangeSearch.h_MatchesCount[0]);



				if ((rangeSearch.h_MatchesCount[0] > 0) && (rangeSearch.h_atomACount[0] > 0))
				{
					gpuErrchk(cudaMemcpy(rangeSearch.h_atomAMatches, rangeSearch.d_atomAMatches, sizeof(int)*rangeSearch.h_MatchesCount[0], cudaMemcpyDeviceToHost));
					gpuErrchk(cudaMemcpy(rangeSearch.h_atomBMatches, rangeSearch.d_atomBMatches, sizeof(int)*rangeSearch.h_MatchesCount[0], cudaMemcpyDeviceToHost));
					int valuesOffset = heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 3 * currentEntry;

					for (int j = 0; j < rangeSearch.h_MatchesCount[0]; j++)
					{
						if (outputType == 2)
						{
							std::cout << "Atom A: " << heldProteinSets.ProteinDataHolder[i].namesSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + rangeSearch.h_atomAMatches[j]] << "\t - Pos : " << rangeSearch.h_atomAMatches[j] << "\t X: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + rangeSearch.h_atomAMatches[j]])) / 1000) << "\t Y: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize + rangeSearch.h_atomAMatches[j]])) / 1000) << "\t Z: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 2 + rangeSearch.h_atomAMatches[j]])) / 1000) << std::endl;
							std::cout << "Atom B: " << heldProteinSets.ProteinDataHolder[i].namesSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + rangeSearch.h_atomBMatches[j]] << "\t - Pos : " << rangeSearch.h_atomBMatches[j] << "\t X: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + rangeSearch.h_atomBMatches[j]])) / 1000) << "\t Y: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize + rangeSearch.h_atomBMatches[j]])) / 1000) << "\t Z: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 2 + rangeSearch.h_atomBMatches[j]])) / 1000) << std::endl << std::endl;

							std::cout << std::endl;
						}
						else if (outputType == 4)
						{
							filePrinter.printLineToOutputFile("Atom A: ", heldProteinSets.ProteinDataHolder[i].namesSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + rangeSearch.h_atomAMatches[j]], "\t - Pos : ", rangeSearch.h_atomAMatches[j], "\t X: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + rangeSearch.h_atomAMatches[j]])) / 1000), "\t Y: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize + rangeSearch.h_atomAMatches[j]])) / 1000), "\t Z: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 2 + rangeSearch.h_atomAMatches[j]])) / 1000));
							filePrinter.printLineToOutputFile("Atom B: ", heldProteinSets.ProteinDataHolder[i].namesSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + rangeSearch.h_atomBMatches[j]], "\t - Pos : ", rangeSearch.h_atomBMatches[j], "\t X: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + rangeSearch.h_atomBMatches[j]])) / 1000), "\t Y: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize + rangeSearch.h_atomBMatches[j]])) / 1000), "\t Z: ", ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 2 + rangeSearch.h_atomBMatches[j]])) / 1000));
							filePrinter.printLineToOutputFile("");
						}

					}


				}
				else
				{
				}

			}


		}

	}

	for (int i = 0; i < 5; i++) //For each of the 5 range lengths of stored protein: Unpin memory.
	{
		proteinsInSet = heldProteinSets.ProteinDataHolder[i].heldEntries;

		if (proteinsInSet > 0)
		{
			gpuErrchk(cudaHostUnregister(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets));
			gpuErrchk(cudaHostUnregister(heldProteinSets.ProteinDataHolder[i].kdTrees));
			gpuErrchk(cudaHostUnregister(heldProteinSets.ProteinDataHolder[i].namesSets));
		}
	}

	gpuErrchk(cudaFree(rangeSearch.d_atomACount));
	gpuErrchk(cudaFree(rangeSearch.d_atomAPositionList));
	gpuErrchk(cudaFree(rangeSearch.d_atomACurrentSearchDimensions));
	gpuErrchk(cudaFree(rangeSearch.d_atomACurrentSearchKdTreePositions));
	gpuErrchk(cudaFree(rangeSearch.d_atomAMatches));
	gpuErrchk(cudaFree(rangeSearch.d_atomBMatches));
	gpuErrchk(cudaFree(rangeSearch.d_MatchesCount));
	gpuErrchk(cudaFree(rangeSearch.d_nextSearchCount));
	gpuErrchk(cudaFree(rangeSearch.d_completionFlag));
	gpuErrchk(cudaFree(rangeSearch.d_kdTreeSets));
	gpuErrchk(cudaFree(rangeSearch.d_namesSets));
	gpuErrchk(cudaFree(rangeSearch.d_xyzCoordsSets));

}





//this could use some work
void calculateInitialBlocksAndThreads(int &blocks, int&threads, int maxEntrySize)
{

	if (maxEntrySize > 512)
	{
		threads = 512;
		blocks = maxEntrySize / threads;
	}
}



__global__ void device_side_ProcessCurrentTreePositionsV3(int* d_xyzValues, short* d_Names, int*  d_kdSetArray, int* d_ViableElementAPairs, int* d_ViableElementBPairs, int* d_ViableElementPairCount, int* d_currentSearchLocations, int*d_currentSearchAElementPos, int* d_currentSearchDimensions, int* CurrentSearchCount, int SizeOfDimensionArray, int maxDistanceSquared, int elementB, int*d_nextSearchCount, int entryNum, int sizeOfKdTree, int maxDistance, int entryBeingWorkedOn)
{
	int id = threadIdx.x + blockDim.x * blockIdx.x;
	/*Variables associated with pointing to correct dataset in large storage array*/
	///////////////////////////
	int currentSetOffset = SizeOfDimensionArray * 3 * entryNum;//As the xyzSet array holds the entries for every single xyz set being processed, the starting position of the current set needs to be added to the equation
	int kdTreeOffSet = sizeOfKdTree*entryNum;
	///////////////////////////


	int threadNamePos = d_currentSearchAElementPos[id]; //This provides the kd position of the A element that all compared items in the search will be measured to.
	int threadkdPos = d_currentSearchLocations[id]; //The position in the kdArray currently being pointed to -- does not include the kdTreeSet offset
	int currentDim = d_currentSearchDimensions[id]; //The dimension of the array, used to choose which dimension to sort by. Queued child processes will recieve currentDim+1

	int currentOffset = SizeOfDimensionArray; //Used for hoping from one dimension in the xyz array to the next -- does not include the xyzSet offset
	int currentLeftChildNode = 0; //Left child node of the current kdTree node -- not yet set;
	int currentRightChildNode = 0; //Left child node of the current kdTree node  -- not yet set ;
	double currentSquaredDistance = 0; //Distance between current node and the sought after node  -- not yet set.
	bool leftChildValid = false; //used to check if left child is worth pursuing
	bool rightChildValid = false; //used to check if right child is worth pursuing
	int ValidChildsOnThisStep = 0; //Maths allows the removal of if statements.
	int resultLocation = id; // This variable acts as a local place holder of the write destination of found matches
	int zOffset = SizeOfDimensionArray * 2; // when switching between the start location of the x,y and z arrays in the global array, y is at SizeOfDimensionArray and z at SizeOfDimensionArray * 2
	int nextSearchSavePosition;

	float elementA_xyz[3]; //Fixed values of the focused element A
	elementA_xyz[0] = d_xyzValues[currentSetOffset + threadNamePos];
	elementA_xyz[1] = d_xyzValues[currentSetOffset + threadNamePos + SizeOfDimensionArray];
	elementA_xyz[2] = d_xyzValues[currentSetOffset + threadNamePos + zOffset];


	float threadxyz[3]; //Holder for the child nodes in the kd tree being looked at. Cuts down on memory calls.
	int currentKdTreePos;

	__syncthreads();


	for (int r = 0; r < 17; r++) //limit is to ensure no more then 14 resulting searches are produced - though no tree is likely to have 25 levels
	{
		currentKdTreePos = d_kdSetArray[kdTreeOffSet + threadkdPos];
		if ((threadkdPos < sizeOfKdTree) && (currentKdTreePos > -1) && (threadkdPos > -1) && (id < CurrentSearchCount[0])) //checking that we aren't looking further then the tree actually extends. 
		{

			threadxyz[0] = d_xyzValues[currentSetOffset + currentKdTreePos];  //Loading the current kd position's coords into the thread
			threadxyz[1] = d_xyzValues[currentSetOffset + currentKdTreePos + SizeOfDimensionArray];
			threadxyz[2] = d_xyzValues[currentSetOffset + currentKdTreePos + zOffset];


			currentLeftChildNode = threadkdPos * 2 + 1;
			currentRightChildNode = threadkdPos * 2 + 2;

			currentSquaredDistance = pow((elementA_xyz[0] - threadxyz[0]), 2) + pow((elementA_xyz[1] - threadxyz[1]), 2) + pow((elementA_xyz[2] - threadxyz[2]), 2); //3 dimensional distance between current and target point

			if (currentSquaredDistance <= maxDistanceSquared)//This section executes if the atoms are close enough for a potential match
			{
				if (d_Names[SizeOfDimensionArray*entryNum + currentKdTreePos] == elementB)//checking if the atom is the sought partner in the pair.
				{
					resultLocation = atomicAdd(&d_ViableElementPairCount[0], 1);
					d_ViableElementAPairs[resultLocation] = threadNamePos;
					d_ViableElementBPairs[resultLocation] = currentKdTreePos;
				}




				if (currentLeftChildNode < sizeOfKdTree)//If atoms are close enough for a match to happen, both child nodes automatically need to be explored - if there are child nodes
				{
					if (currentRightChildNode < sizeOfKdTree)
					{
						nextSearchSavePosition = atomicAdd(&d_nextSearchCount[0], 1);

						d_currentSearchLocations[nextSearchSavePosition] = currentRightChildNode;
						d_currentSearchAElementPos[nextSearchSavePosition] = threadNamePos;
						d_currentSearchDimensions[nextSearchSavePosition] = (currentDim + 1) % 3;
					}
					threadkdPos = currentLeftChildNode;
					currentDim = (currentDim + 1) % 3;
				}
				else if (currentRightChildNode < sizeOfKdTree)
				{
					threadkdPos = currentRightChildNode;
					currentDim = (currentDim + 1) % 3;
				}
				else //If both child nodes return invalid, we have reached the end of the search of this branch.
				{
					return;
				}
			}
			else
			{
				ValidChildsOnThisStep = 0;
				leftChildValid = 0;
				rightChildValid = 0;
				int targetElementInCurrentDim = elementA_xyz[currentDim];
				int currentThreadxyz = threadxyz[currentDim];
				int absoluteDistanceBetweenCurrentAndTargetNode = abs(currentThreadxyz - targetElementInCurrentDim);
				bool currentNodeLessThenTargetInCurrentDim = (currentThreadxyz < targetElementInCurrentDim);
				bool  currentNodeGreaterThenTargetInCurrentDim = (currentThreadxyz > targetElementInCurrentDim);

				int d_kdSetArrayLeftChildNode = d_kdSetArray[kdTreeOffSet + currentLeftChildNode];
				int currentNodeOutsideRequiredRangeInGivenDimension = (absoluteDistanceBetweenCurrentAndTargetNode > maxDistance);
				int currentD_xyzValues;
				if (d_kdSetArrayLeftChildNode > -1)
				{
					currentD_xyzValues = d_xyzValues[currentSetOffset + d_kdSetArrayLeftChildNode + currentOffset*currentDim];
					int currentLeftChildNodeFurtherFromTargetThenCurrentNodeInGivenDimension = (abs(currentD_xyzValues - targetElementInCurrentDim) > absoluteDistanceBetweenCurrentAndTargetNode);
					int currentLeftNodeInSameDirectionAsCurrentNodeFromTargetNode = ((currentNodeLessThenTargetInCurrentDim && (currentD_xyzValues<targetElementInCurrentDim)) || (currentNodeGreaterThenTargetInCurrentDim && (currentD_xyzValues>targetElementInCurrentDim)));
					int AllConditionsMetForLeftNodeInvalid = currentNodeOutsideRequiredRangeInGivenDimension*currentLeftChildNodeFurtherFromTargetThenCurrentNodeInGivenDimension*currentLeftNodeInSameDirectionAsCurrentNodeFromTargetNode; //These can be put straight into the booleans if they confirmed to work;
					leftChildValid = 1 - AllConditionsMetForLeftNodeInvalid;
				}

				int d_kdSetArrayRightChildNode = d_kdSetArray[kdTreeOffSet + currentRightChildNode];
				if (d_kdSetArrayRightChildNode > -1)
				{
					currentD_xyzValues = d_xyzValues[currentSetOffset + d_kdSetArrayRightChildNode + currentOffset*currentDim];
					int currentRightChildNodeFurtherFromTargetThenCurrentNodeInGivenDimension = (abs(currentD_xyzValues - targetElementInCurrentDim) > absoluteDistanceBetweenCurrentAndTargetNode);
					int currentRightNodeInSameDirectionAsCurrentNodeFromTargetNode = ((currentNodeLessThenTargetInCurrentDim && (currentD_xyzValues<targetElementInCurrentDim)) || (currentNodeGreaterThenTargetInCurrentDim && (currentD_xyzValues>targetElementInCurrentDim)));
					int AllConditionsMetForRightNodeInvalid = currentNodeOutsideRequiredRangeInGivenDimension*currentRightChildNodeFurtherFromTargetThenCurrentNodeInGivenDimension*currentRightNodeInSameDirectionAsCurrentNodeFromTargetNode; //These can be put straight into the booleans if they confirmed to work;
					rightChildValid = 1 - AllConditionsMetForRightNodeInvalid;
				}


				ValidChildsOnThisStep = 0 + leftChildValid + rightChildValid;

				if (ValidChildsOnThisStep == 0)
				{
					return; //End Point for thread.
				}
				else if (ValidChildsOnThisStep == 1)
				{
					currentDim = (currentDim + 1) % 3;

					if (leftChildValid == true)
					{
						threadkdPos = currentLeftChildNode;
					}
					else
					{
						threadkdPos = currentRightChildNode;
					}
				}
				else
				{
					nextSearchSavePosition = atomicAdd(&d_nextSearchCount[0], 1);
					d_currentSearchLocations[nextSearchSavePosition] = currentRightChildNode;
					d_currentSearchAElementPos[nextSearchSavePosition] = threadNamePos;
					d_currentSearchDimensions[nextSearchSavePosition] = (currentDim + 1) % 3;

					threadkdPos = currentLeftChildNode;
					currentDim = (currentDim + 1) % 3;
				}
			}

		}
		else
		{
			return;
		}
	}




	return;
}


__global__ void SetNextCountAsCurrentAndCheckFlag(int* d_currentSearchCount, int* d_nextSearchCount, int* d_completionFlag)
{
	if (d_nextSearchCount[0] == 0)
	{
		d_completionFlag[0] = 1;
	}
	else
	{
		d_currentSearchCount[0] = d_nextSearchCount[0];
		d_nextSearchCount[0] = 0;
	}
	return;
}
