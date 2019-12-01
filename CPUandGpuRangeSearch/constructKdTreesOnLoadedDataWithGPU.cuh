#ifndef GPUKDTREERANGESEARCHWITH
#define GPUKDTREERANGESEARCHWITH
#include "cuda_runtime.h"


#include "printToFileHandler.h"
#include "miscFunctions.h"
#include <string>
#include <vector>

#include <stdlib.h>
#include <iostream>
#include <stdio.h>
//#include "dataStructures.h"
#include "ESBTL_ProteinHandler.h"
#include "SearchStructureConstruction.h"

void constructKdTreesOnLoadedDataOnCPUWithGPUSorting(ProteinDataHandler &ProteinData);

void Parallel_Device_Side_bitonic_sort(int *h_coordArray1, int * h_reverseValuePositionArray, int* d_forwardPositionArray, int*d_reversePositionArray, int * d_coordsHolder, int *h_coordArray2, int * h_reverseValuePositionArray2, int* d_forwardPositionArray2, int*d_reversePositionArray2, int * d_coordsHolder2, int *h_coordArray3, int * h_reverseValuePositionArray3, int* d_forwardPositionArray3, int*d_reversePositionArray3, int * d_coordsHolder3, int numOfBlock, int numOfThreads);
__global__ void resetReverseAndForwardArrays(int* d_xbackArray, int* d_xforwardArray, int* d_ybackArray, int* d_yforwardArray, int* d_zbackArray, int* d_zforwardArray, int length);

void kdTree_constructorLocal(int* kdArray, int* xPositionReferenceList, int* xBackPositionReferenceList, int* yPositionReferenceList, int* yBackPositionReferenceList, int*  zPositionReferenceList, int*  zBackPositionReferenceList, int  numOfRealInputElements, int lengthOfHolderArrays);

void constructKdTreesOnLoadedDataOnGPU(ProteinDataHandler &ProteinData);
void Device_Side_kdTree_constructor(int* d_kdArray, int* d_xPositionReferenceList, int* d_xBackPositionReferenceList, int* d_yPositionReferenceList, int* d_yBackPositionReferenceList, int*  d_zPositionReferenceList, int*  d_zBackPositionReferenceList, int  *numOfRealInputElements, int lengthOfHolderArrays);

__global__ void online_kdConstructor_master_node(int * d_kdArray, int * d_xPositionReferenceList, int * d_xBackPositionReferenceList, int * d_yPositionReferenceList, int * d_yBackPositionReferenceList, int * d_zPositionReferenceList, int * d_zBackPositionReferenceList, int numberOfCoords);
__global__ void online_kdConstructor_slave_nodeX(int * d_xBackPositionReferenceList, int totalNumOfConsecutiveTasks, int* kdArray, int* d_xPositionReferenceList, int*d_yPositionReferenceList, int* d_zPositionReferenceList);
__global__ void online_kdConstructor_slave_nodeY(int * d_yBackPositionReferenceList, int totalNumOfConsecutiveTasks, int* kdArray, int* d_xPositionReferenceList, int*d_yPositionReferenceList, int* d_zPositionReferenceList);
__global__ void online_kdConstructor_slave_nodeZ(int * d_zBackPositionReferenceList, int totalNumOfConsecutiveTasks, int* kdArray, int* d_xPositionReferenceList, int*d_yPositionReferenceList, int* d_zPositionReferenceList);

void GPUTreeWithCPULoad(ProteinDataHandler &ProteinData);
#endif
