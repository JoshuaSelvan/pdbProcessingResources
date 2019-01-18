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

		d_array[i]=d_array[i]+5;
	

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




//for small sets
void gpuKdTreeUnoptimisedRangeSearchAllLoadedSets(rangeSearchSettings &settings, ProteinDataHandler heldProteinSets, AtomToNumHashTable atomReferenceTable)
{
	//all data has been preloaded into the host side ProteinDataHandler object. All that is needed is to initialise containers, moves sets of data to the gpu, process those sets of data and then return the results.

	std::cout << "PERFORMING GPU BASED KD-TREE RANGE SEARCH" << std::endl;
	//Initialise host and device data holders
	gpuRangeSearchResources rangeSearch;
	int safeHolderSize = 32000 *140* sizeof(int)*3;
	
	int blocks = 1;
	int threads = 1;
	int currentEntry = 0;
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
		cudaHostRegister(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets, proteinsInSet*entrySize*sizeof(int) * 3, 0);
		cudaHostRegister(heldProteinSets.ProteinDataHolder[i].namesSets, proteinsInSet*entrySize *sizeof(short), 0);
		cudaHostRegister(heldProteinSets.ProteinDataHolder[i].kdTrees, proteinsInSet*kdEntrySize*sizeof(int), 0);
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
			if (heldProteinSets.ProteinDataHolder[i].heldEntries>0)
		{
std::cout<< proteinsInSet*entrySize*sizeof(int) * 3<<std::endl;	
	gpuErrchk(cudaMemcpy(rangeSearch.d_xyzCoordsSets, heldProteinSets.ProteinDataHolder[i].xyzCoordsSets, proteinsInSet*entrySize*sizeof(int) * 3, cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpyAsync(rangeSearch.d_namesSets, heldProteinSets.ProteinDataHolder[i].namesSets, proteinsInSet*entrySize *sizeof(short), cudaMemcpyHostToDevice/*, streams[0]*/));
		gpuErrchk(cudaMemcpyAsync(rangeSearch.d_kdTreeSets, heldProteinSets.ProteinDataHolder[i].kdTrees, proteinsInSet*kdEntrySize*sizeof(int), cudaMemcpyHostToDevice/*, streams[0]*/));

		std::cout << "Processing set: " << i << " on the gpu" << std::endl;
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
			blocks = heldProteinSets.ProteinDataHolder[i].MaxEntrySize/512;
			device_side__locateElement << <blocks, threads/*, 0, streams[0] */>> >(rangeSearch.d_namesSets, currentEntry, soughtAtomANumber, rangeSearch.d_atomAPositionList, rangeSearch.d_atomACount, rangeSearch.d_atomACurrentSearchDimensions, rangeSearch.d_atomACurrentSearchKdTreePositions, heldProteinSets.ProteinDataHolder[i].MaxEntrySize);

			gpuErrchk(cudaMemcpyAsync(rangeSearch.h_atomACount, rangeSearch.d_atomACount, sizeof(int), cudaMemcpyDeviceToHost, streams[0]));
			//gpuErrchk(cudaMemcpyAsync(rangeSearch.h_atomAPositionList, rangeSearch.d_atomAPositionList, safeHolderSize, cudaMemcpyDeviceToHost, streams[1]));
			std::cout << "File: " << currentEntry << " - atom A count: " << rangeSearch.h_atomACount[0] << "  ";
			
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
					device_side_ProcessCurrentTreePositionsV3 << <rangeSearch.blocks, rangeSearch.threads >> >(rangeSearch.d_xyzCoordsSets, rangeSearch.d_namesSets, rangeSearch.d_kdTreeSets, rangeSearch.d_atomAMatches, rangeSearch.d_atomBMatches, rangeSearch.d_MatchesCount, rangeSearch.d_atomACurrentSearchKdTreePositions, rangeSearch.d_atomAPositionList, rangeSearch.d_atomACurrentSearchDimensions, rangeSearch.d_atomACount, heldProteinSets.ProteinDataHolder[i].MaxEntrySize, maxDistanceSquared, soughtAtomBNumber, rangeSearch.d_nextSearchCount, currentEntry, heldProteinSets.ProteinDataHolder[i].KdTreeSize, settings.requiredProximity, currentEntry);





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
			if ((rangeSearch.h_MatchesCount[0] > 0) && (rangeSearch.h_atomACount[0] > 0))
			{
				
				std::cout << "Matches in file: "<< rangeSearch.h_MatchesCount[0] << std::endl;
				gpuErrchk(cudaMemcpy(rangeSearch.h_atomAMatches, rangeSearch.d_atomAMatches, sizeof(int)*rangeSearch.h_MatchesCount[0], cudaMemcpyDeviceToHost));
				gpuErrchk(cudaMemcpy(rangeSearch.h_atomBMatches, rangeSearch.d_atomBMatches, sizeof(int)*rangeSearch.h_MatchesCount[0], cudaMemcpyDeviceToHost));
				//std::cout << std::endl << "Matches in file " << currentEntry << std::endl;
				//std::cout << "Number of matches found: " << rangeSearch.h_MatchesCount[0] << std::endl << "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_" << std::endl;
				int valuesOffset = heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 3 * currentEntry;
				for (int j = 0; j < rangeSearch.h_MatchesCount[0]; j++)
				{
				//	std::cout << "Atom A: " << heldProteinSets.ProteinDataHolder[i].namesSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + rangeSearch.h_atomAMatches[j]] << "\t - Pos : " << rangeSearch.h_atomAMatches[j] << "\t X: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + rangeSearch.h_atomAMatches[j]])) / 1000) << "\t Y: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize + rangeSearch.h_atomAMatches[j]])) / 1000) << "\t Z: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 2 + rangeSearch.h_atomAMatches[j]])) / 1000) << std::endl;
				//	std::cout << "Atom B: " << heldProteinSets.ProteinDataHolder[i].namesSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + rangeSearch.h_atomBMatches[j]] << "\t - Pos : " << rangeSearch.h_atomBMatches[j] << "\t X: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + rangeSearch.h_atomBMatches[j]])) / 1000) << "\t Y: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize + rangeSearch.h_atomBMatches[j]])) / 1000) << "\t Z: " << ((double(heldProteinSets.ProteinDataHolder[i].xyzCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 2 + rangeSearch.h_atomBMatches[j]])) / 1000) << std::endl << std::endl;
		//temp outtake		//	std::cout << "Atom A: " << heldProteinSets.ProteinDataHolder[i].namesSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + rangeSearch.h_atomAMatches[j]] << "\t - Pos : " << rangeSearch.h_atomAMatches[j] << "\t X: " << ((double(heldProteinSets.ProteinDataHolder[i].xCoordsSets[valuesOffset + rangeSearch.h_atomAMatches[j]])) / 1000) << "\t Y: " << ((double(heldProteinSets.ProteinDataHolder[i].yCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize + rangeSearch.h_atomAMatches[j]])) / 1000) << "\t Z: " << ((double(heldProteinSets.ProteinDataHolder[i].zCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 2 + rangeSearch.h_atomAMatches[j]])) / 1000) << std::endl;
				//	std::cout << "Atom B: " << heldProteinSets.ProteinDataHolder[i].namesSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + rangeSearch.h_atomBMatches[j]] << "\t - Pos : " << rangeSearch.h_atomBMatches[j] << "\t X: " << ((double(heldProteinSets.ProteinDataHolder[i].xCoordsSets[valuesOffset + rangeSearch.h_atomBMatches[j]])) / 1000) << "\t Y: " << ((double(heldProteinSets.ProteinDataHolder[i].yCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize + rangeSearch.h_atomBMatches[j]])) / 1000) << "\t Z: " << ((double(heldProteinSets.ProteinDataHolder[i].zCoordsSets[valuesOffset + heldProteinSets.ProteinDataHolder[i].MaxEntrySize * 2 + rangeSearch.h_atomBMatches[j]])) / 1000) << std::endl << std::endl;

				}
				//std::cout << "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_" << std::endl;
			}
			else
				std::cout << "Matches in file: 0" << std::endl;
		}
		}

	}
	//Need to finish moving stuff here.

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


	for (int r = 0; r < 25; r++) //limit is to ensure no more then 14 resulting searches are produced - though no tree is likely to have 25 levels
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
