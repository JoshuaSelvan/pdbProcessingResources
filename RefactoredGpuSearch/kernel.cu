
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include "miscFunctions.h"
#include "dataStructures.h"
#include "rangeSearch.h"
#include <time.h>
#include <string>


int main()
{
	clock_t start, stop;
	rangeSearchSettings currentSettings;
	std::string runSettingsFile = "runDetails.txt";

	start = clock();
	loadRunDetails(runSettingsFile, currentSettings);

	initiateRangeSearch(currentSettings);
	stop = clock();


	print_elapsed(start, stop, "Total run time:");
	return 0;

}
