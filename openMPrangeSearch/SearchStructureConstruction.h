#ifndef SEARCHSTRUCTURECONSTRUCTION
#define SEARCHSTRUCTURECONSTRUCTION

#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include "dataStructures.h"
#include "ESBTL_ProteinHandler.h"


//handles the construction and use of generic search structures.

void kdTree_constructor(int* kdArray, int* xPositionReferenceList, int* xBackPositionReferenceList, int* yPositionReferenceList, int* yBackPositionReferenceList, int*  zPositionReferenceList, int*  zBackPositionReferenceList, int  numOfRealInputElements, int lengthOfHolderArrays);

void populate_X_KdNode(int taskNo, kdTask * pendingList, int * d_xBackPositionReferenceList, int totalNumOfConsecutiveTasks, int* kdArray, int* d_xPositionReferenceList, int*d_yPositionReferenceList, int* d_zPositionReferenceList);
void populate_Y_KdNode(int taskNo, kdTask * pendingList, int * d_yBackPositionReferenceList, int totalNumOfConsecutiveTasks, int* kdArray, int* d_xPositionReferenceList, int*d_yPositionReferenceList, int* d_zPositionReferenceList);
void populate_Z_KdNode(int taskNo, kdTask * pendingList, int * d_zBackPositionReferenceList, int totalNumOfConsecutiveTasks, int* kdArray, int* d_xPositionReferenceList, int*d_yPositionReferenceList, int* d_zPositionReferenceList);


void sortReferenceArrayByQuickSort(int *tempValueArray, int* BackCoordsLoader, int startPos, int endPos);
void setforwardReferenceArrayFromBack(int *ForwardPosLoader, int *xBackCoordsLoader, int numOfAtoms);

void formatSecondaryPositionStructure(short* namesSetsArray, int maximumLengthOfChains, int EntryCount, int* chainLengthsList, int* compositionCountsList, int*compositionListArray, int*compositionPointerArray);
void searchEntryInSecondaryPositionStructureForAtom(short soughtAtomNum, int EntryNum, int maxSizeOfChains, int* compositionCountsList, int*compositionListArray, int*compositionPointerArray, int*outPutArray, int*outputCount);

int sizeOfFinalKdLevel(int chainSize, int &numberOfLevels);


void constructKdTreesOnLoadedDataOnCPU(ProteinDataHandler &ProteinData);



void processSingleFile(rangeSearchSettings &settings, ProteinDataHandler heldProteinSets, AtomToNumHashTable atomReferenceTable, int rangeToProcess, int entryToProcess);
void constructKdTreesOnLoadedDataOnCPUOMP(ProteinDataHandler &ProteinData); 

//void setSecondarySearchStructure(ProteinDataHandler &ProteinData);

#endif
