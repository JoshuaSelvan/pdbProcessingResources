
#include "miscFunctions.h"
#include "dataStructures.h"



void switchLoadingAndProcessingSets(int &loading, int &processing)
{
	int temp = loading;
	loading = processing;
	processing = temp;
	return;
}





void print_chrono_elapsed(std::chrono::time_point<std::chrono::high_resolution_clock> start, std::chrono::time_point<std::chrono::high_resolution_clock> end, std::string caption)
{
    typedef std::chrono::high_resolution_clock Time;
    typedef std::chrono::milliseconds ms;
    typedef std::chrono::duration<float> fsec;

    fsec fs = end-start;

    ms d = std::chrono::duration_cast<ms>(fs);
       std::cout<<caption<<" " << d.count() << "ms\n";
    }
    
     


void print3ArraysToFile(int* arrayA, int* arrayB, int* arrayC, int length)
{
	std::ofstream sourceFileLoader;
	sourceFileLoader.open("SavedCoordinates.txt");
	
	for (int i = 0; i < length; i++)
	{
		sourceFileLoader << arrayA[i] << "\t" << arrayB[i] << "\t" << arrayC[i] << "\n";

	}
	sourceFileLoader.close();


}



void print1ArraysToFile(int* arrayA, int length)
{
	std::ofstream sourceFileLoader;
	sourceFileLoader.open("SavedCoordinates.txt");

	for (int i = 0; i < length; i++)
	{
		sourceFileLoader<<i<<"\t" << arrayA[i] << "\n";

	}
	sourceFileLoader.close();


}


void resetMultiRunDetails()
{

	std::ofstream sourceFileLoader;
	sourceFileLoader.open("multiRunTimer.txt");
	std::string holderString;

	sourceFileLoader << "Number_of_timing_sets_to_process:" << std::endl;
	sourceFileLoader << "2" << std::endl;
	sourceFileLoader << "Set1:" << std::endl;
	sourceFileLoader << "Search_Type:(1_cpuBruteForce,2_cpuKdTree,3_gpuBruteForce,4_gpuBruteForce)" << std::endl;
	sourceFileLoader << "0" << std::endl;
	sourceFileLoader << "atoms:" << std::endl;
	sourceFileLoader << "C1'" << std::endl;
	sourceFileLoader << "C2'" << std::endl;
	sourceFileLoader << "range_to_Search:" << std::endl;
	sourceFileLoader << "0" << std::endl;
	sourceFileLoader << "number_of_proximities_to_test:" << std::endl;
	sourceFileLoader << "4" << std::endl;
	sourceFileLoader << "proximities_list:" << std::endl;
	sourceFileLoader << "0.5" << std::endl;
	sourceFileLoader << "1" << std::endl;
	sourceFileLoader << "4" << std::endl;
	sourceFileLoader << "8" << std::endl;
	sourceFileLoader << "number_of_entry_ranges:" << std::endl;
	sourceFileLoader << "4" << std::endl;
	sourceFileLoader << "entry_ranges:" << std::endl;
	sourceFileLoader << "200" << std::endl;
	sourceFileLoader << "400" << std::endl;
	sourceFileLoader << "800" << std::endl;
	sourceFileLoader << "1600" << std::endl;
	sourceFileLoader << "-------------------" << std::endl;
	sourceFileLoader << "Set2:" << std::endl;
	sourceFileLoader << "atoms:" << std::endl;
	sourceFileLoader << "CB" << std::endl;
	sourceFileLoader << "CG" << std::endl;
	sourceFileLoader << "range_to_Search:" << std::endl;
	sourceFileLoader << "0" << std::endl;
	sourceFileLoader << "number_of_proximities_to_test:" << std::endl;
	sourceFileLoader << "4" << std::endl;
	sourceFileLoader << "proximities_list:" << std::endl;
	sourceFileLoader << "0.5" << std::endl;
	sourceFileLoader << "1" << std::endl;
	sourceFileLoader << "4" << std::endl;
	sourceFileLoader << "8" << std::endl;
	sourceFileLoader << "number_of_entry_ranges:" << std::endl;
	sourceFileLoader << "4" << std::endl;
	sourceFileLoader << "entry_ranges:" << std::endl;
	sourceFileLoader << "200" << std::endl;
	sourceFileLoader << "400" << std::endl;
	sourceFileLoader << "800" << std::endl;
	sourceFileLoader << "1600" << std::endl;
	sourceFileLoader.close();


}

int loadMultiRunDetails(std::vector<multiRunDetailsSet>& multiSet, std::string fileName)
{
	std::ifstream sourceFileLoader;
	sourceFileLoader.open(fileName.c_str());
	std::string holderString;
	int numberOfEntries;
	int numberOfProximitesToTest;
	float nextProximity;
	int numberOfEntriesToSearch;
	int nextnumberOfEntries;

	sourceFileLoader >> holderString;
	sourceFileLoader >> numberOfEntries;
	sourceFileLoader >> holderString;

	for (int i = 0; i < numberOfEntries; i++)
	{
		
		multiRunDetailsSet currentset;
		currentset.proximitiesCount = 0;
		sourceFileLoader >> holderString;
		sourceFileLoader >> currentset.searchType;
		sourceFileLoader >> holderString;
		sourceFileLoader >> currentset.atom1; 
		sourceFileLoader >> currentset.atom2;
		sourceFileLoader >> holderString;
		sourceFileLoader >> currentset.currentRange;
		sourceFileLoader >> holderString;
		sourceFileLoader >> numberOfProximitesToTest;
		sourceFileLoader >> holderString;

		for (int i = 0; i < numberOfProximitesToTest; i++)
		{
			sourceFileLoader >> nextProximity;
			currentset.proximities.push_back(nextProximity);
			currentset.proximitiesCount++;
		}
		sourceFileLoader >> holderString;
		sourceFileLoader >> numberOfEntriesToSearch;
		sourceFileLoader >> holderString;
		for (int i = 0; i < numberOfEntriesToSearch; i++)
		{
			sourceFileLoader >> nextnumberOfEntries;
			currentset.entriesToSearchEachRound.push_back(nextnumberOfEntries);
		}

		sourceFileLoader >> holderString;
		sourceFileLoader >> holderString;

		multiSet.push_back(currentset);
	}
	
	sourceFileLoader.close();
	return numberOfEntries;
}


void loadRunDetails(std::string sourceFile, rangeSearchSettings &currentSettings)
{

std::ifstream sourceFileLoader;
sourceFileLoader.open(sourceFile.c_str());
std::string holderString;


sourceFileLoader >> holderString;
sourceFileLoader>>holderString;
sourceFileLoader>>holderString;
sourceFileLoader>>currentSettings.inputListFileLocation;
sourceFileLoader >> holderString;
sourceFileLoader >> currentSettings.searchType;
sourceFileLoader >> holderString;
sourceFileLoader >> currentSettings.debugLevel;
sourceFileLoader >> holderString;
sourceFileLoader >> currentSettings.resultsPrintFormat;


sourceFileLoader>>holderString;
sourceFileLoader>>currentSettings.numberOfFiles;
sourceFileLoader>>holderString;
sourceFileLoader>>currentSettings.AtomTypeOne;
sourceFileLoader>>holderString;
sourceFileLoader>>currentSettings.AtomTypeTwo;
sourceFileLoader>>holderString;
sourceFileLoader>>currentSettings.requiredProximity;
currentSettings.requiredProximity = currentSettings.requiredProximity*1000;
sourceFileLoader>>holderString;    
sourceFileLoader >> currentSettings.maxProteinLengthPerRange[0];
sourceFileLoader>>holderString;
sourceFileLoader>>currentSettings.maxNumOfEntriesPerRange[0];
sourceFileLoader>>holderString;
sourceFileLoader>>currentSettings.maxProteinLengthPerRange[1];
sourceFileLoader>>holderString;
sourceFileLoader>>currentSettings.maxNumOfEntriesPerRange[1];
sourceFileLoader>>holderString;
sourceFileLoader>>currentSettings.maxProteinLengthPerRange[2];
sourceFileLoader>>holderString;
sourceFileLoader>>currentSettings.maxNumOfEntriesPerRange[2];
sourceFileLoader>>holderString;
sourceFileLoader>>currentSettings.maxProteinLengthPerRange[3];
sourceFileLoader>>holderString;
sourceFileLoader>>currentSettings.maxNumOfEntriesPerRange[3];
sourceFileLoader>>holderString;
sourceFileLoader>>currentSettings.maxProteinLengthPerRange[4];
sourceFileLoader>>holderString;
sourceFileLoader>>currentSettings.maxNumOfEntriesPerRange[4];

sourceFileLoader.close();    

};



void print_elapsed(clock_t start, clock_t end, std::string captionBeforeTimeDifference)
{
    
	double elapsed = ((double)(end - start)) / CLOCKS_PER_SEC;
	std::cout<<captionBeforeTimeDifference<<" "<<elapsed<<std::endl;
     
}

 

//This function creates a default input file for people to put their settings into.
void resetInputTextFile()
{
	std::ofstream sourceFileLoader;
	sourceFileLoader.open("runDetails.txt");
	std::string holderString;


	sourceFileLoader << "CPU_range_search_parameters"<<std::endl;
	sourceFileLoader << "---------------------------" << std::endl;
	sourceFileLoader << "Protein_file_list_location:" << std::endl;
	sourceFileLoader << "C:\\Users\\Joshua\\Desktop\\refactoredThesisExecutalbes\\RefactoredRangeSearchV1\\RefactoredRangeSearchV1\\pdbFileSet\\pdbListShort.txt" << std::endl;
	sourceFileLoader << "Search_Type:(1-bruteForce,2-kdTree,3-gpuBruteForce,4-gpuKdTree)" << std::endl;
	sourceFileLoader << "1" << std::endl;
	sourceFileLoader << "Number_of_files:" << std::endl;
	sourceFileLoader << "10000" << std::endl;
	sourceFileLoader << "Atom_A:" << std::endl;
	sourceFileLoader << "C1" << std::endl;
	sourceFileLoader << "Atom_B:" << std::endl;
	sourceFileLoader << "C2" << std::endl;
	sourceFileLoader << "Required_Proximity:" << std::endl;
	sourceFileLoader << "4" << std::endl;
	sourceFileLoader << "Range_1_max_protein_length:" << std::endl;
	sourceFileLoader << "1024" << std::endl;
	sourceFileLoader << "Range_1_max_entry_count:" << std::endl;
	sourceFileLoader << "1000" << std::endl;
	sourceFileLoader << "Range_2_max_protein_length:" << std::endl;
	sourceFileLoader << "2048" << std::endl;
	sourceFileLoader << "Range_2_max_entry_count:" << std::endl;
	sourceFileLoader << "2500" << std::endl;
	sourceFileLoader << "Range_3_max_protein_length:" << std::endl;
	sourceFileLoader << "4096" << std::endl;
	sourceFileLoader << "Range_3_max_entry_count:" << std::endl;
	sourceFileLoader << "750" << std::endl;
	sourceFileLoader << "Range_4_max_protein_length:" << std::endl;
	sourceFileLoader << "8192" << std::endl;
	sourceFileLoader << "Range_4_max_entry_count:" << std::endl;
	sourceFileLoader << "200" << std::endl;
	sourceFileLoader << "Range_5_max_protein_length:" << std::endl;
	sourceFileLoader << "16384" << std::endl;
	sourceFileLoader << "Range_5_max_entry_count:" << std::endl;
	sourceFileLoader << "10" << std::endl;

	sourceFileLoader.close();


}


void reformatFileList(std::string unformattedFile, std::string newPrefix, std::string newSuffix)
{
	std::ifstream sourceFileLoader;
	sourceFileLoader.open(unformattedFile.c_str());

	std::string holderString;

	std::ofstream destFile;
	destFile.open("reformattedFileList.txt");

	while (sourceFileLoader >> holderString)
	{
		destFile << newPrefix << holderString << newSuffix << std::endl;
	}

	sourceFileLoader.close();
	destFile.close();

}


int largestSizeInArrayPos(int *theArray, int arraySize)
{
	int largest = 0;
	for (int i = 0; i < arraySize; i++)
	{
		if (theArray[i]>theArray[largest])
			largest = i;
	}

	return largest;
}


void loadMultiBatchRangeSizes(std::string sourceFile, rangeSearchSettings &currentSettings)
{

	std::ifstream sourceFileLoader;
	sourceFileLoader.open(sourceFile.c_str());
	std::string holderString;


	sourceFileLoader >> holderString;
	sourceFileLoader >> holderString;
	
	sourceFileLoader >> currentSettings.multiBatchRangeSizes[0];
	sourceFileLoader >> holderString;
	sourceFileLoader >> currentSettings.multiBatchRangeSizes[1];
	sourceFileLoader >> holderString;
	sourceFileLoader >> currentSettings.multiBatchRangeSizes[2];
	sourceFileLoader >> holderString;
	sourceFileLoader >> currentSettings.multiBatchRangeSizes[3];
	sourceFileLoader >> holderString;
	sourceFileLoader >> currentSettings.multiBatchRangeSizes[4];
	
	sourceFileLoader.close();

};



void resetMultiBatchRangeSizes()
{
	std::ofstream sourceFileLoader;
	sourceFileLoader.open("MultiBatchSizes.txt");
	std::string holderString;

	sourceFileLoader << "Set_Ranges_For_Multi_Batch_GPU_Searches:" << std::endl;
	sourceFileLoader << "Range_0:" << std::endl;
	sourceFileLoader << "2" << std::endl;
	sourceFileLoader << "Range_1:" << std::endl;
	sourceFileLoader << "2" << std::endl;
	sourceFileLoader << "Range_2:" << std::endl;
	sourceFileLoader << "2" << std::endl;
	sourceFileLoader << "Range_3:" << std::endl;
	sourceFileLoader << "2" << std::endl;
	sourceFileLoader << "Range_4:" << std::endl;
	sourceFileLoader << "1" << std::endl;
	
	sourceFileLoader.close();

}
