#include <mpi.h>

#include <stdlib.h>
#include <iostream>
#include <string>
#include <cmath>
#include <fstream>

#include <stdio.h>
#include "miscFunctions.h"
#include "dataStructures.h"
#include "rangeSearchHandler.h"
#include <time.h>

#include <utmpx.h>


using namespace std;






int main(int argc, char** argv) {

	// Initialize the MPI environment
	MPI_Init(NULL, NULL);
	//Find out rank, size
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);


	char processor_name[MPI_MAX_PROCESSOR_NAME];
    	int name_len;
    	MPI_Get_processor_name(processor_name, &name_len);



	if(world_rank==0)
		cout<<"World size is: "<<world_size<<endl;

	
	rangeSearchSettings currentSettings;
	std::string runSettingsFile = "runDetails.txt";
	loadRunDetails(runSettingsFile, currentSettings);
	string inputLocation;
	float totalInputFiles =currentSettings.numberOfFiles;
	int personalStartingPointInList=0;




	if(world_rank==0)
		cout<<"Total number of input files is: "<<totalInputFiles<<endl;



	int numberOfFilesPerProcessor = (int) ceil( (float)totalInputFiles / world_size ); 
	int personalFilesToProcess=0;
	if(world_rank<(world_size-1))
	{
		personalFilesToProcess= (int) ceil( (float)totalInputFiles / world_size );
		personalStartingPointInList=numberOfFilesPerProcessor*world_rank;
	}
	else
	{
		personalFilesToProcess= (int) floor( (float)totalInputFiles / world_size );
		personalStartingPointInList=numberOfFilesPerProcessor*world_rank;
	}


	char numstr[21];


	string outputFileName = "ResultFromRank";
	sprintf(numstr, "%d", world_rank);
	outputFileName += numstr;
	outputFileName += ".txt";


	//printf( "cpu = %d\n", sched_getcpu() );
	cout<<"ProcessorRank: "<<world_rank<<" on host: "<<processor_name <<" Is processing: "<<personalFilesToProcess<<" starting at position: "<<personalStartingPointInList<<" and is printing to: "<<outputFileName<<endl;


	clock_t  start,end;
	start =clock();
	initiateRangeSearch(currentSettings,personalFilesToProcess,personalStartingPointInList);


	end=clock();
	print_elapsed(start, end,"total run time: ");


	MPI_Finalize();
return 0;
}








void quickSort(int* arr, int left, int right) {
	int i = left, j = right;
	int tmp;
	int pivot = arr[(left + right) / 2];

	/* partition */
	while (i <= j) {
		while (arr[i] < pivot)
			i++;
		while (arr[j] > pivot)
			j--;
		if (i <= j) {
			tmp = arr[i];
			arr[i] = arr[j];
			arr[j] = tmp;
			i++;
			j--;
		}
	};

	/* recursion */
	if (left < j)
		quickSort(arr, left, j);
	if (i < right)
		quickSort(arr, i, right);
}

void extendedQuickSort(int* arrayOfValues,int*backwardReference, int left, int right) {
	int i = left, j = right;
	int tmp;
	int tmp2;
	int pivot = arrayOfValues[(left + right) / 2];

	/* partition */
	while (i <= j) {
		while (arrayOfValues[i] < pivot)
			i++;
		while (arrayOfValues[j] > pivot)
			j--;
		if (i <= j) {

			tmp = arrayOfValues[i];
			tmp2 = backwardReference[i];

			arrayOfValues[i] = arrayOfValues[j];
			backwardReference[i] = backwardReference[j];


			arrayOfValues[j] = tmp;
			backwardReference[j] = tmp2;
			i++;
			j--;
		}
	};

	/* recursion */
	if (left < j)
		extendedQuickSort(arrayOfValues, backwardReference, left, j);
	if (i < right)
		extendedQuickSort(arrayOfValues,  backwardReference, i, right);
}

