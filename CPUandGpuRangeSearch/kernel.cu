
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include "miscFunctions.h"
#include "dataStructures.h"
#include "rangeSearch.h"
#include <time.h>
#include <string>
#include "gpuMemoryLimitTester.cuh"

int main()
{

	//checkGpuMemoryLimit();
	//return;


		clock_t start, stop;
	rangeSearchSettings currentSettings;
	std::string runSettingsFile = "runDetails.txt";

	//reformatFileList("C:\\Users\\Joshua\\Desktop\\refactoredThesisExecutalbes\\RefactoredGpuSearch\\RefactoredGpuSearch\\pdbFileSet\\pdbListMedium.txt", "pdbFileSet\\", ""); - not up to date. Do not use!
	//resetMultiRunDetails();//working




		start = clock();
	loadRunDetails(runSettingsFile, currentSettings);
	initiateRangeSearch(currentSettings);

	//ProteinDataHandler heldProteinSets(currentSettings);
	//heldProteinSets.createSubListsContainingValidFilesForEachRangeOnly(currentSettings);

	
		stop = clock();


	print_elapsed(start, stop, "Run time including load time: ");
	return 0;

}
