
#include "miscFunctions.h"
#include "dataStructures.h"

void loadRunDetails(std::string sourceFile, rangeSearchSettings &currentSettings)
{

std::ifstream sourceFileLoader;
sourceFileLoader.open(sourceFile/*.c_str()*/);
std::string holderString;


sourceFileLoader >> holderString;
sourceFileLoader>>holderString;
sourceFileLoader>>holderString;
sourceFileLoader>>currentSettings.inputListFileLocation;
sourceFileLoader >> holderString;
sourceFileLoader >> currentSettings.searchType;
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
	std::cout<<captionBeforeTimeDifference<<" "<<elapsed;
     
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