#include "SearchStructureConstruction.h"



void populate_X_KdNode(int taskNo, kdTask * pendingList, int * d_xBackPositionReferenceList, int totalNumOfConsecutiveTasks, /*int searchRange,*/ int* kdArray, int* d_xPositionReferenceList, int*d_yPositionReferenceList, int* d_zPositionReferenceList)
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

void populate_Y_KdNode(int taskNo, kdTask * pendingList, int * d_yBackPositionReferenceList, int totalNumOfConsecutiveTasks, /*int searchRange,*/ int* kdArray, int* d_xPositionReferenceList, int*d_yPositionReferenceList, int* d_zPositionReferenceList)
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

void populate_Z_KdNode(int taskNo, kdTask * pendingList, int * d_zBackPositionReferenceList, int totalNumOfConsecutiveTasks, /*int searchRange,*/ int* kdArray, int* d_xPositionReferenceList, int*d_yPositionReferenceList, int* d_zPositionReferenceList)
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


void searchEntryInSecondaryPositionStructureForAtom(short soughtAtomNum, int EntryNum, int maxSizeOfChains, int* compositionCountsList, int*compositionListArray, int*compositionPointerArray, int*outPutArray, int*outputCount)
{
	outputCount[0] = 0;
	int startPosOffSet = EntryNum*maxSizeOfChains;
	int inEntryAtomTypeNum = -1;
	int resultsPositionOffSet = 0;
	int numOfMatches = 0;
	for (int i = 0; i < compositionCountsList[EntryNum]; i++)
	{
		if (compositionListArray[startPosOffSet + i*2] == soughtAtomNum)
		{
			numOfMatches = compositionListArray[startPosOffSet + i*2 + 1];
			outputCount[0] = numOfMatches;
			for (int j = 0; j < numOfMatches; j++)
			{
				outPutArray[j] = compositionPointerArray[startPosOffSet + resultsPositionOffSet + j];
			}
			i = i + compositionCountsList[EntryNum] + 1;
		}
		else
		{
			resultsPositionOffSet = resultsPositionOffSet + compositionListArray[startPosOffSet + i*2 + 1];
			if (compositionListArray[startPosOffSet + i] == -1)
				std::cout << "Error occured" << std::endl;
		}

	}


}

int sizeOfFinalKdLevel(int chainSize, int &numberOfLevels)
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

void sortReferenceArrayByQuickSort(int* tempValueArray, int* BackCoordsLoader, int startPos, int endPos)
{
	int i = startPos, j = endPos;
	int tmp;
	int tmp2;
	int pivot = tempValueArray[(startPos + endPos) / 2];

	/* partition */
	while (i <= j) {
		while (tempValueArray[i] < pivot)
			i++;
		while (tempValueArray[j] > pivot)
			j--;
		if (i <= j) {

			tmp = tempValueArray[i];
			tmp2 = BackCoordsLoader[i];

			tempValueArray[i] = tempValueArray[j];
			BackCoordsLoader[i] = BackCoordsLoader[j];


			tempValueArray[j] = tmp;
			BackCoordsLoader[j] = tmp2;
			i++;
			j--;
		}
	};

	/* recursion */
	if (startPos < j)
		sortReferenceArrayByQuickSort(tempValueArray, BackCoordsLoader, startPos, j);
	if (i < endPos)
		sortReferenceArrayByQuickSort(tempValueArray, BackCoordsLoader, i, endPos);

}

void setforwardReferenceArrayFromBack(int *ForwardPosLoader, int *xBackCoordsLoader, int numOfAtoms)
{
	for (int i = 0; i < numOfAtoms; i++)
	{
		ForwardPosLoader[xBackCoordsLoader[i]] = i;
	}
}

void kdTree_constructor(int* kdArray, int* xPositionReferenceList, int* xBackPositionReferenceList, int* yPositionReferenceList, int* yBackPositionReferenceList, int*  zPositionReferenceList, int*  zBackPositionReferenceList, int  numOfRealInputElements, int lengthOfHolderArrays)
{
	int numOfLevels;
	sizeOfFinalKdLevel(numOfRealInputElements, numOfLevels);
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
				populate_X_KdNode(j, taskList, xBackPositionReferenceList, activeTasks, kdArray, xPositionReferenceList, yPositionReferenceList, zPositionReferenceList);
		else if (i % 3 == 1)
			for (int j = 0; j < activeTasks; j++)
				populate_Y_KdNode(j, taskList, yBackPositionReferenceList, activeTasks, kdArray, xPositionReferenceList, yPositionReferenceList, zPositionReferenceList);
		else if (i % 3 == 2)
			for (int j = 0; j < activeTasks; j++)
				populate_Z_KdNode(j, taskList, zBackPositionReferenceList, activeTasks, kdArray, xPositionReferenceList, yPositionReferenceList, zPositionReferenceList);
		activeTasks = activeTasks * 2;
	}
}




void constructKdTreesOnLoadedDataOnCPU(ProteinDataHandler &ProteinData) 
{
	dataLoader kdTreeConstructionArtifacts;


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



	for (int i = 0; i<5; i++)
	{

		for (int j = 0; j<ProteinData.ProteinDataHolder[i].heldEntries; j++)
		{	
			if(j%500==0)
				std::cout<<"-";
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
			std::fill(kdTreeConstructionArtifacts.kdTreeholder, kdTreeConstructionArtifacts.kdTreeholder + ProteinData.ProteinDataHolder[i].MaxEntrySize*2, -1);



			sortReferenceArrayByQuickSort(kdTreeConstructionArtifacts.xValues, kdTreeConstructionArtifacts.xBackPositionReferenceList, 0, (ProteinData.ProteinDataHolder[i].proteinLengthCounts[j] - 1));
			setforwardReferenceArrayFromBack(kdTreeConstructionArtifacts.xForwardPositionReferenceList, kdTreeConstructionArtifacts.xBackPositionReferenceList, ProteinData.ProteinDataHolder[i].proteinLengthCounts[j]);

			sortReferenceArrayByQuickSort(kdTreeConstructionArtifacts.yValues, kdTreeConstructionArtifacts.yBackPositionReferenceList, 0, (ProteinData.ProteinDataHolder[i].proteinLengthCounts[j] - 1));
			setforwardReferenceArrayFromBack(kdTreeConstructionArtifacts.yForwardPositionReferenceList, kdTreeConstructionArtifacts.yBackPositionReferenceList, ProteinData.ProteinDataHolder[i].proteinLengthCounts[j]);

			sortReferenceArrayByQuickSort(kdTreeConstructionArtifacts.zValues, kdTreeConstructionArtifacts.zBackPositionReferenceList, 0, (ProteinData.ProteinDataHolder[i].proteinLengthCounts[j] - 1));
			setforwardReferenceArrayFromBack(kdTreeConstructionArtifacts.zForwardPositionReferenceList, kdTreeConstructionArtifacts.zBackPositionReferenceList, ProteinData.ProteinDataHolder[i].proteinLengthCounts[j]);



			//construct kdTree
			kdTree_constructor(kdTreeConstructionArtifacts.kdTreeholder, kdTreeConstructionArtifacts.xForwardPositionReferenceList, kdTreeConstructionArtifacts.xBackPositionReferenceList, kdTreeConstructionArtifacts.yForwardPositionReferenceList, kdTreeConstructionArtifacts.yBackPositionReferenceList, kdTreeConstructionArtifacts.zForwardPositionReferenceList, kdTreeConstructionArtifacts.zBackPositionReferenceList, ProteinData.ProteinDataHolder[i].proteinLengthCounts[j], ProteinData.ProteinDataHolder[i].MaxEntrySize);


			//move kdTree to main storage arrays
			for (int y = 0; y<ProteinData.ProteinDataHolder[i].MaxEntrySize * 2; y++)
			{
				ProteinData.ProteinDataHolder[i].kdTrees[y + j* ProteinData.ProteinDataHolder[i].KdTreeSize] = kdTreeConstructionArtifacts.kdTreeholder[y];

			}
		}
	}
}


