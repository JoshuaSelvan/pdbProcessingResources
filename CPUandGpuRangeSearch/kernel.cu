
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include "miscFunctions.h"
#include "dataStructures.h"
#include "rangeSearch.h"
#include <time.h>
#include <string>
#include "gpuMemoryLimitTester.cuh"
#include "interface.h"
#include <iostream>



int main()
{

	displayProgramDetails();
	int proceed = checkUserOperations();

	if (proceed == 0)
		return 0;
	//checkGpuMemoryLimit();
	//return;


	clock_t start, stop;
	rangeSearchSettings currentSettings;
	std::string runSettingsFile = "runDetails.txt";
	
	//This line is redundant. Type in 4 at begin to select the rename option
	//reformatFileList("range1FileList.txt","/local/pdb/",".gz", "internalRange1list.txt");
	//return 0;	


	//resetMultiRunDetails();//working




		start = clock();
	loadRunDetails(runSettingsFile, currentSettings);
	initiateRangeSearch(currentSettings);

	//ProteinDataHandler heldProteinSets(currentSettings);
	//heldProteinSets.createSubListsContainingValidFilesForEachRangeOnly(currentSettings);

	
		stop = clock();


	print_elapsed(start, stop, "Total run time: ");
	return 0;

}
