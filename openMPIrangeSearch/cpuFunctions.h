#ifndef GPUKDTREEFUNCTIONS
#define GPUKDTREEFUNCTIONS

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
using namespace std;



struct kdTask {
	int kdDestination;
	int minX;
	int maxX;
	int minY;
	int maxY;
	int minZ;
	int maxZ;
};



void print_elapsed(clock_t start, clock_t stop)
{
	double elapsed = ((double)(stop - start)) / CLOCKS_PER_SEC;
	printf("Elapsed time: %.3fs\n", elapsed);
}

void setIntArrayValues(int *array, int value, int lengthOfArray)
{
	
		for (int j = 0; j < lengthOfArray; j++)
		{
			array[j] = value;
		}
	
};
void setShortArrayValues(short *array, int value, int lengthOfArray)
{

	for (int j = 0; j < lengthOfArray; j++)
	{
		array[j] = value;
	}

};
void printStateOfSingleIntArray(int* array, int numOfAtoms, int MaxEntrySize)
{
	ofstream arrayChecker;
	arrayChecker.open("BitonicInputs.txt", ios::app);
	arrayChecker << "numOfAtoms: " << numOfAtoms << endl << "MaxEntrySize: " << MaxEntrySize << endl;

	for (int j = 0; j < MaxEntrySize; j++)
	{
		arrayChecker << (j) << "\t\t" << array[j] << endl;
	}
	arrayChecker << "//////////////////////////////////////////////////////////////////////////////////////" << endl << "//////////////////////////////////////////////////////////////////////////////////////" << endl << "//////////////////////////////////////////////////////////////////////////////////////" << endl;
	arrayChecker.close();
}
void printStateOfDeviceSideKdTreesHolder(int* kdTreesArray, int kdTreeSize, int numOfEntries)//function cant work from the file, mearly a blueprint;
{

	ofstream arrayChecker;
	arrayChecker.open("kdTreesPrintOut.txt", ios::app);
	arrayChecker << "kdTreeSize: " << kdTreeSize << endl << "numOfEntries: " << numOfEntries << endl;

	int* h_kdTreesArray = (int*)malloc(sizeof(int)*kdTreeSize*numOfEntries);

	for (int j = 0; j < kdTreeSize*numOfEntries; j++)
	{
		arrayChecker << (j) << "\t\t" << (j) % kdTreeSize << h_kdTreesArray[j] << endl;
		if ((j + 1) % kdTreeSize == 0)
			arrayChecker << "//////////////////////////////////////////////////////////////////////////////////////" << endl << "//////////////////////////////////////////////////////////////////////////////////////" << endl << "//////////////////////////////////////////////////////////////////////////////////////" << endl;

	}
	arrayChecker.close();

}
void printStateOfInputsToBitonicSort(int *h_coordArray, /*int* numOfAtomsInEntry, */int * h_reverseValuePositionArray, int numOfBlock, int numOfThreads, int* d_forwardPositionArray, int*d_reversePositionArray, int * d_coordsHolder) //This function cant actually work so im just storing the blueprint here
{
	ofstream arrayChecker;
	arrayChecker.open("BitonicInputs.txt", ios::app);
	arrayChecker << "numOfBlock: " << numOfBlock << endl << "numOfThreads: " << numOfThreads << endl;
	int * tempA = (int*)malloc(numOfBlock*numOfThreads*sizeof(int));
	int * tempB = (int*)malloc(numOfBlock*numOfThreads*sizeof(int));
	int * tempC = (int*)malloc(numOfBlock*numOfThreads*sizeof(int));

	for (int j = 0; j < numOfBlock*numOfThreads; j++)
	{
		arrayChecker << (j) << "\t\t" << h_coordArray[j] << "\t\t" << h_reverseValuePositionArray << "\t\t" << tempA << "\t\t" << tempB<<"\t\t"<<tempC<<endl;
	}
	arrayChecker << "//////////////////////////////////////////////////////////////////////////////////////" << endl << "//////////////////////////////////////////////////////////////////////////////////////" << endl << "//////////////////////////////////////////////////////////////////////////////////////" << endl;
	arrayChecker.close();
}
void printStateOfAllHostXYZCoordsToFile(int *xCoords, int*yCoords, int*zCoords, int Entries, int MaxEntrySize)
{
	ofstream arrayChecker;
	arrayChecker.open("stateOfInput.txt", ios::app);
	arrayChecker << "Entrys: " <<Entries<<endl<<"MaxEntrySize: "<<MaxEntrySize<< endl;

		
		for (int j = 0; j < Entries*MaxEntrySize; j++)
		{
			arrayChecker << (j) << "\t\t" << j%MaxEntrySize << "\t\t" << ((double(xCoords[j])) / 1000) << "\t\t" << ((double(yCoords[j])) / 1000) << "\t\t" << ((double(zCoords[j]))/1000) << endl;
			if ((j+1)%MaxEntrySize==0)
				arrayChecker << "//////////////////////////////////////////////////////////////////////////////////////" << endl << "//////////////////////////////////////////////////////////////////////////////////////" << endl << "//////////////////////////////////////////////////////////////////////////////////////" << endl;

		}
		arrayChecker.close();
}
void printStateOfHostXYZandXandYandZToFile(int *xCoords, int*yCoords, int*zCoords, int* xyzCoords, int Entries, int MaxEntrySize)
{
	ofstream arrayChecker;
	arrayChecker.open("stateOfInput.txt", ios::app);
	arrayChecker << "Entrys: " << Entries << endl << "MaxEntrySize: " << MaxEntrySize << endl;
	int p = 0;
	
	for (int j = 0; j < Entries*MaxEntrySize; j++)
	{
		arrayChecker << (j) << "\t\t" << j%MaxEntrySize << "\t\t" << ((double(xCoords[j])) / 1000) << "\t\t" << ((double(xyzCoords[p])) / 1000) << "\t\t" << ((double(yCoords[j])) / 1000) << "\t\t" << ((double(xyzCoords[p + MaxEntrySize])) / 1000) << "\t\t" << ((double(zCoords[j])) / 1000) << "\t\t" << ((double(xyzCoords[p + MaxEntrySize*2])) / 1000) << endl;
		
		p++;
		if ((j + 1) % MaxEntrySize == 0)
		{
			p = p + MaxEntrySize*2;
			arrayChecker << "//////////////////////////////////////////////////////////////////////////////////////" << endl << "//////////////////////////////////////////////////////////////////////////////////////" << endl << "//////////////////////////////////////////////////////////////////////////////////////" << endl;
		}
	}
	arrayChecker.close();
}

void printStateOfSingleShortArray(short* array, int numOfAtoms, int MaxEntrySize)
{
	ofstream arrayChecker;
	arrayChecker.open("BitonicInputs.txt", ios::app);
	arrayChecker << "numOfAtoms: " << numOfAtoms << endl << "MaxEntrySize: " << MaxEntrySize << endl;

	for (int j = 0; j < MaxEntrySize; j++)
	{
		arrayChecker << (j) << "\t\t" << array[j] << endl;
	}
	arrayChecker << "//////////////////////////////////////////////////////////////////////////////////////" << endl << "//////////////////////////////////////////////////////////////////////////////////////" << endl << "//////////////////////////////////////////////////////////////////////////////////////" << endl;
	arrayChecker.close();
}

void printStateOfCurrentBackAndForthReferenceTables(int *backRef, int*forwardRef, int arraySize, int entrycount)
{
	ofstream arrayChecker;
	arrayChecker.open("stateOfInput.txt", ios::app);
	arrayChecker << "arraySize: " << arraySize << endl << "entrycount: " << entrycount << endl;

	for (int j = 0; j < arraySize; j++)
	{
		arrayChecker << (j) << "\t\t" << backRef[j] << "\t\t" << forwardRef[j] << endl;
	}
	arrayChecker << "//////////////////////////////////////////////////////////////////////////////////////" << endl << "//////////////////////////////////////////////////////////////////////////////////////" << endl << "//////////////////////////////////////////////////////////////////////////////////////" << endl;
	arrayChecker.close();
}


void printStateOfReferenceAndReferencedArrays(int *backRef, int*xCoords, int*yCoords, int*zCoords, short*names, int arraySize, int entryNum,int numOfAtoms)
{
	ofstream arrayChecker;
	arrayChecker.open("stateOfInput.txt", ios::app);
	arrayChecker << "arraySize: " << arraySize << endl << "entryNum: " << entryNum << endl<<"NumOfAtoms: "<<numOfAtoms<<endl;


	for (int j = 0; j < arraySize; j++)
	{
		arrayChecker << (j) << "\t\t" << backRef[j] << "\t\t" << "\t\t" << ((double(xCoords[backRef[j]])) / 1000) << "\t\t" << ((double(yCoords[backRef[j]])) / 1000) << "\t\t" << ((double(zCoords[backRef[j]])) / 1000) <<"\t\t"<< names[backRef[j]] << endl;
	}
	arrayChecker << "//////////////////////////////////////////////////////////////////////////////////////" << endl << "//////////////////////////////////////////////////////////////////////////////////////" << endl << "//////////////////////////////////////////////////////////////////////////////////////" << endl;
	arrayChecker.close();
}

void setChainLengthAttributesFromFile(string source, int& lengthOfChain, int& numOfRealInputElements, int& threads, int& blocks)
{
	ifstream inputfile(source.c_str());
	inputfile >> numOfRealInputElements;


	lengthOfChain = 2;
	int counter = 1;
	while (numOfRealInputElements > lengthOfChain)
	{
		lengthOfChain = lengthOfChain * 2;
		counter++;
	}

	if (lengthOfChain > 512)
	{
		threads = 512;
		blocks = lengthOfChain / threads;
	}
	else
	{
		if (counter % 2 == 0)
		{
			threads = int(pow(2, ((double)counter / 2)));
			blocks = lengthOfChain / threads;
		}
		else
		{
			threads = pow(2, ((double)(counter + 1) / 2));
			blocks = lengthOfChain / threads;
		}
	}
	inputfile.close();

}
void ImportTextFileData(string source, int lengthOfChain, int& numOfRealInputElements, vector<string>& Names, float * xvals, float * yvals, float * zvals)
{
	ifstream inputfile(source.c_str());
	string inputStringHolder;
	float inputStreamHolder;

	inputfile >> numOfRealInputElements;
	for (int i = 0; i < numOfRealInputElements; i++)
	{
		inputfile >> inputStreamHolder;
		xvals[i] = inputStreamHolder;
	}
	for (int i = 0; i < numOfRealInputElements; i++)
	{

		inputfile >> inputStreamHolder;
		yvals[i] = inputStreamHolder;
	}
	for (int i = 0; i < numOfRealInputElements; i++)
	{
		inputfile >> inputStreamHolder;
		zvals[i] = inputStreamHolder;
	}
	for (int i = 0; i < numOfRealInputElements; i++)
	{
		inputfile >> inputStringHolder; Names.push_back(inputStringHolder);
	}

	for (int i = numOfRealInputElements; i < lengthOfChain; i++)
	{
		xvals[i] = 9999;
		yvals[i] = 9999;
		zvals[i] = 9999;
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






int sizeOfInitialGrid(int chainSize)
{
	int gridSize = 2;
	if (chainSize > gridSize)
		gridSize = gridSize * 2;

	return gridSize;


}






void offline_kdConstructor_slave_nodeX(int taskNo, kdTask * pendingList, int * d_xBackPositionReferenceList, int totalNumOfConsecutiveTasks, /*int searchRange,*/ int* kdArray, int* d_xPositionReferenceList, int*d_yPositionReferenceList, int* d_zPositionReferenceList)
{



	kdTask task = pendingList[taskNo];

	if (task.kdDestination < 0)
		return;

	int activeElementCount = 0;
	int secondCount = -1;
	int elementToBeExtracted = -1;
	int posOfExtractTarget = -1;

	for (int i = task.minX; i < task.maxX; i++)
	{
		if (d_xPositionReferenceList[d_xBackPositionReferenceList[i]] != -1)
		{
			if (d_yPositionReferenceList[d_xBackPositionReferenceList[i]] <= task.maxY && d_yPositionReferenceList[d_xBackPositionReferenceList[i]] >= task.minY)
			{
				if (d_zPositionReferenceList[d_xBackPositionReferenceList[i]] <= task.maxZ && d_zPositionReferenceList[d_xBackPositionReferenceList[i]] >= task.minZ)
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


	if (activeElementCount % 2 == 0)
		elementToBeExtracted = activeElementCount / 2 - 1;
	else
		elementToBeExtracted = (activeElementCount - 1) / 2;

	int i = task.minX;
	while (posOfExtractTarget == -1)
	{
		if (d_xPositionReferenceList[d_xBackPositionReferenceList[i]] != -1)
			if (d_yPositionReferenceList[d_xBackPositionReferenceList[i]]<=task.maxY && d_yPositionReferenceList[d_xBackPositionReferenceList[i]] >= task.minY)
				if (d_zPositionReferenceList[d_xBackPositionReferenceList[i]]<=task.maxZ && d_zPositionReferenceList[d_xBackPositionReferenceList[i]] >= task.minZ)
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

}
void offline_kdConstructor_slave_nodeY(int taskNo, kdTask * pendingList, int * d_yBackPositionReferenceList, int totalNumOfConsecutiveTasks, /*int searchRange,*/ int* kdArray, int* d_xPositionReferenceList, int*d_yPositionReferenceList, int* d_zPositionReferenceList)
{

	kdTask task = pendingList[taskNo];

	int activeElementCount = 0;
	int secondCount = -1;
	int elementToBeExtracted = -1;
	int posOfExtractTarget = -1;
	for (int i = task.minY; i < task.maxY; i++) 
	{
		if (d_yPositionReferenceList[d_yBackPositionReferenceList[i]] != -1)
			if (d_xPositionReferenceList[d_yBackPositionReferenceList[i]]<=task.maxX && d_xPositionReferenceList[d_yBackPositionReferenceList[i]] >= task.minX)
				if (d_zPositionReferenceList[d_yBackPositionReferenceList[i]]<=task.maxZ && d_zPositionReferenceList[d_yBackPositionReferenceList[i]] >= task.minZ)
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
		if (d_yPositionReferenceList[d_yBackPositionReferenceList[i]] != -1)
			if (d_xPositionReferenceList[d_yBackPositionReferenceList[i ]]<=task.maxX && d_xPositionReferenceList[d_yBackPositionReferenceList[i ]] >= task.minX)
				if (d_zPositionReferenceList[d_yBackPositionReferenceList[i]]<=task.maxZ && d_zPositionReferenceList[d_yBackPositionReferenceList[i ]] >= task.minZ)
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
void offline_kdConstructor_slave_nodeZ(int taskNo, kdTask * pendingList, int * d_zBackPositionReferenceList, int totalNumOfConsecutiveTasks, /*int searchRange,*/ int* kdArray, int* d_xPositionReferenceList, int*d_yPositionReferenceList, int* d_zPositionReferenceList)
{

	kdTask task = pendingList[taskNo];
	
	if (task.kdDestination < 0)
		return;

	int activeElementCount = 0;
	int secondCount = -1;
	int elementToBeExtracted = -1;
	int posOfExtractTarget = -1;

	for (int i = task.minZ; i < task.maxZ; i++) 
	{
		if (d_zPositionReferenceList[d_zBackPositionReferenceList[i]] != -1)
		{
			if (d_yPositionReferenceList[d_zBackPositionReferenceList[i]] <= task.maxY && d_yPositionReferenceList[d_zBackPositionReferenceList[i]] >= task.minY)
			{
				if (d_xPositionReferenceList[d_zBackPositionReferenceList[i]] <= task.maxX && d_xPositionReferenceList[d_zBackPositionReferenceList[i]] >= task.minX)
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
			if (d_zPositionReferenceList[d_zBackPositionReferenceList[i]] != -1)
				if (d_yPositionReferenceList[d_zBackPositionReferenceList[i]] <= task.maxY && d_yPositionReferenceList[d_zBackPositionReferenceList[i]] >= task.minY)
					if (d_xPositionReferenceList[d_zBackPositionReferenceList[i]] <= task.maxX && d_xPositionReferenceList[d_zBackPositionReferenceList[i]] >= task.minX)
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
void offline_kdTree_constructor(int* kdArray, int* xPositionReferenceList, int* xBackPositionReferenceList, int* yPositionReferenceList, int* yBackPositionReferenceList, int*  zPositionReferenceList, int*  zBackPositionReferenceList, int  numOfRealInputElements, int lengthOfHolderArrays)
{
	int numOfLevels;
	int taskArraySize = sizeOfFinalKdLevel(numOfRealInputElements, numOfLevels); 
	kdTask *taskList = (kdTask*)malloc(sizeof(kdTask)*(taskArraySize*2));
	

	int numberOfCoords = numOfRealInputElements;
	for (int i = 0; i < taskArraySize * 2;i++)
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



	for (int i = 0; i < (numOfLevels); i++)
	{
	
		if (i % 3 == 0)
			for (int j = 0; j < activeTasks; j++)
				offline_kdConstructor_slave_nodeX(j, taskList, xBackPositionReferenceList, activeTasks, /*numberOfCoords,*/ kdArray, xPositionReferenceList, yPositionReferenceList, zPositionReferenceList);
		else if (i % 3 == 1)
			for (int j = 0; j < activeTasks; j++)
				offline_kdConstructor_slave_nodeY(j, taskList, yBackPositionReferenceList, activeTasks, /*numberOfCoords,*/ kdArray, xPositionReferenceList, yPositionReferenceList, zPositionReferenceList);
		else if (i % 3 == 2)
			for (int j = 0; j < activeTasks; j++)
				offline_kdConstructor_slave_nodeZ(j, taskList, zBackPositionReferenceList, activeTasks, /*numberOfCoords,*/ kdArray, xPositionReferenceList, yPositionReferenceList, zPositionReferenceList);
		activeTasks = activeTasks * 2;
	}


}




float vocalDistance3D(int atomAX, int atomAY, int atomAZ, float x, float y, float z)
{
	float temp = pow((atomAX - x), 2) + pow((atomAY - y), 2) + pow((atomAZ - z), 2);
	cout << pow((atomAX - x), 2) << endl;
	cout << pow((atomAY - y), 2) << endl;
	cout << pow((atomAZ - z), 2) << endl;
	
	return temp;
}

int distance3D(int atomAX, int atomAY, int atomAZ, int atomBx, int atomBy, float atomBz)
{
	int temp = pow((atomAX - atomBx), 2) + pow((atomAY - atomBy), 2) + pow((atomAZ - atomBz), 2);
	return temp;
}

bool withinGraceZone(int atomAX, int atomAY, int atomAZ, int xCurrent, int yCurrent, int zCurrent, int maxDistance)
{
	if ((pow((atomAX - xCurrent), 2) + pow((atomAY - yCurrent), 2) + pow((atomAZ - zCurrent), 2)) < (sqrt(3 * pow(maxDistance, 2))))
		return true;
	else
		return false;

}

bool viableDirection(int atomAX, int atomAY, int atomAZ,  int xCurrent, int yCurrent, int zCurrent, int xNext, int yNext, int zNext, int treeLevel, int requiredDistance)
{
	int dim = treeLevel % 3;

	if (dim == 0)
		if (atomAX > xCurrent)
			if ((xCurrent > xNext)&&(abs(xCurrent-atomAX)>requiredDistance))
				return false;
			else
				return true;
		else
			if ((xCurrent < xNext) && (abs(xCurrent - atomAX)>requiredDistance))
				return false;
			else
				return true;
	else if (dim == 1)
		if (atomAY > yCurrent)
			if ((yCurrent > yNext) && (abs(yCurrent - atomAY)>requiredDistance))
				return false;
			else
				return true;
		else
			if ((yCurrent < yNext) && (abs(yCurrent - atomAY)>requiredDistance))
				return false;
			else
				return true;
	else if (dim == 2)
		if (atomAZ > zCurrent)
			if ((zCurrent > zNext) && (abs(zCurrent - atomAZ)>requiredDistance))
				return false;
			else
				return true;
		else
			if ((zCurrent < zNext) && (abs(zCurrent - atomAZ)>requiredDistance))
				return false;
			else
				return true;
	return false; 
}

float distance1D(int atomAX, int atomAY, int atomAZ, int atomBx, int atomBy, float atomBz, int treeLevel)
{
	int dim = treeLevel % 3;

	if (dim == 0)
		return abs(atomAX - atomBx);
	else if (dim == 1)
		return abs(atomAY - atomBy);
	else if (dim == 2)
		return abs(atomAZ - atomBz);

	return 0;
}

void rangeSearchCpu(int atomAX, int atomAY, int atomAZ, int elementTwo, int maxDistance, vector<int> &resultsVector, short *names, int *xValueArray, int *yValueArray, int*zValueArray, int *kdArray, int kdArraySize, int TreePos, int treeLevel, int* runcount, int valueArraysStartPositionOffset, int kdTreeStartPositionOffset)
{


	
	int currentElement = names[valueArraysStartPositionOffset + kdArray[kdTreeStartPositionOffset + TreePos]];

	runcount[0] = runcount[0]++;
	if (kdArray[kdTreeStartPositionOffset+TreePos] == -1)
		int q = 7;

	else
	{

		int leftChildNode = TreePos * 2 + 1;
		int rightChildNode = TreePos * 2 + 2;
		float distanceFromAtomAToLocation = distance3D(atomAX, atomAY, atomAZ, xValueArray[valueArraysStartPositionOffset + kdArray[kdTreeStartPositionOffset + TreePos]], yValueArray[valueArraysStartPositionOffset + kdArray[kdTreeStartPositionOffset + TreePos]], zValueArray[valueArraysStartPositionOffset + kdArray[kdTreeStartPositionOffset + TreePos]]);


		int q = leftChildNode + 7;
		if (distanceFromAtomAToLocation < pow(maxDistance, 2)) {

			if (names[valueArraysStartPositionOffset + kdArray[kdTreeStartPositionOffset + TreePos]] == (elementTwo))
			{
				
				resultsVector.push_back(kdArray[kdTreeStartPositionOffset+TreePos]);
			}
	
			if (leftChildNode < kdArraySize)
				if (kdArray[kdTreeStartPositionOffset+leftChildNode] != -1)
					rangeSearchCpu(atomAX, atomAY, atomAZ, elementTwo, maxDistance, resultsVector, names, xValueArray, yValueArray, zValueArray, kdArray, kdArraySize, leftChildNode, (treeLevel + 1), runcount, valueArraysStartPositionOffset, kdTreeStartPositionOffset);
			if (rightChildNode < kdArraySize)
				if (kdArray[kdTreeStartPositionOffset+rightChildNode] != -1)
					rangeSearchCpu(atomAX, atomAY, atomAZ, elementTwo, maxDistance, resultsVector, names, xValueArray, yValueArray, zValueArray, kdArray, kdArraySize, rightChildNode, (treeLevel + 1), runcount, valueArraysStartPositionOffset, kdTreeStartPositionOffset);
		}
		else
		{
			if (leftChildNode < kdArraySize)
				if ((viableDirection(atomAX, atomAY, atomAZ, xValueArray[valueArraysStartPositionOffset + kdArray[kdTreeStartPositionOffset + TreePos]], yValueArray[valueArraysStartPositionOffset + kdArray[kdTreeStartPositionOffset + TreePos]], zValueArray[valueArraysStartPositionOffset + kdArray[kdTreeStartPositionOffset + TreePos]], xValueArray[valueArraysStartPositionOffset + kdArray[kdTreeStartPositionOffset + leftChildNode]], yValueArray[valueArraysStartPositionOffset + kdArray[kdTreeStartPositionOffset + leftChildNode]], zValueArray[valueArraysStartPositionOffset + kdArray[kdTreeStartPositionOffset + leftChildNode]], treeLevel,maxDistance) == true)); 
                if (kdArray[kdTreeStartPositionOffset+leftChildNode] != -1)
						rangeSearchCpu(atomAX, atomAY, atomAZ, elementTwo, maxDistance, resultsVector, names, xValueArray, yValueArray, zValueArray, kdArray, kdArraySize, leftChildNode, (treeLevel + 1), runcount, valueArraysStartPositionOffset, kdTreeStartPositionOffset);
	
			if (rightChildNode < kdArraySize)
				if ((viableDirection(atomAX, atomAY, atomAZ, xValueArray[valueArraysStartPositionOffset + kdArray[kdTreeStartPositionOffset + TreePos]], yValueArray[valueArraysStartPositionOffset + kdArray[kdTreeStartPositionOffset + TreePos]], zValueArray[valueArraysStartPositionOffset + kdArray[kdTreeStartPositionOffset + TreePos]], xValueArray[valueArraysStartPositionOffset + kdArray[kdTreeStartPositionOffset + rightChildNode]], yValueArray[valueArraysStartPositionOffset + kdArray[kdTreeStartPositionOffset + rightChildNode]], zValueArray[valueArraysStartPositionOffset + kdArray[kdTreeStartPositionOffset + rightChildNode]], treeLevel, maxDistance) == true))
					if (kdArray[kdTreeStartPositionOffset+rightChildNode] != -1)
						rangeSearchCpu(atomAX, atomAY, atomAZ, elementTwo, maxDistance, resultsVector, names, xValueArray, yValueArray, zValueArray, kdArray, kdArraySize, rightChildNode, (treeLevel + 1), runcount, valueArraysStartPositionOffset, kdTreeStartPositionOffset);
		}
	}

}

void primitiveLocateAtomA(short soughtAtom, short* namesSets, int lengthOfEntries, int EntryNumber, int* atomAPosList, int* atomACount)
{
	atomACount[0] = 0;
	for (int i = 0; i < lengthOfEntries; i++)
	{
		if (namesSets[lengthOfEntries*EntryNumber + i] == soughtAtom)
		{
			atomAPosList[atomACount[0]] = lengthOfEntries*EntryNumber + i;
			atomACount[0]++;

		}


	}


}







#endif

