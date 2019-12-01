#include "interface.h"

void displayProgramDetails()
{
	std::cout << std::endl;
	std::cout << "Salutations!" <<std::endl;
	std::cout << "This program has been written to perform range searches on PDB files utilizing brute force and kd-tree algorithms on either a CPU or a GPU."<<std::endl;
std::cout << "This program accepts a list of the locations of a batch of unzipped ent PDB files and then processes them in accordance to the settings chosen in the runDetails.txt file" << std::endl <<std::endl;

	std::cout << "[1] To recreate the runDetails.txt file."<<std::endl;
	std::cout << "[2] To run the program with the options currently written in RunDetails.txt"<<std::endl;
	std::cout << "[3] To abort"<<std::endl; 
}

int checkUserOperations()
{
	int input;
	std::cin >> input;

	if( input == 1)
	{
		recreateRunDetails();
		return 0;
	}
	else if (input == 2)
	{
		return 1;
	}	
	else if (input == 4)
	{
	    std::string fileName;
 	    std::string pre;
	    std::string suf;
	    std::string neww;
	    std::cout << "Input file to convert"<<std::endl;

	    std::cin >> fileName;

	    std::cout << "Input prefix"<<std::endl;

    	    std::cin >> pre;

            std::cout << "Input sufix"<<std::endl;

            std::cin >> suf;

            std::cout << "Input new file Name"<<std::endl;

            std::cin >> neww;


            reformatFileList(fileName,pre,suf,neww);
            return 0;
	}





else
	{
		return 0;
	}

}


void recreateRunDetails()
{
        std::ofstream sourceFileLoader;
        sourceFileLoader.open("runDetails.txt");
        std::string holderString;

        sourceFileLoader << "Range_search_parameters"<<std::endl;
        sourceFileLoader << "---------------------------" << std::endl;
        sourceFileLoader << "Protein_file_list_location:" << std::endl;
        sourceFileLoader << "pdbListShort.txt" << std::endl;
        sourceFileLoader << "Search_Type:(1-bruteForce,2-gpuBruteForceSimple,3-gpuBruteForceSingleload,4-gpuBruteForceBatchLoad,5-kdtree,6-gpuKdTreeSimple,7-gpuKdTreeSingleLoad,8-gpuKdTreeBatchLoad,9-MultiProcessRun,10-kdTreeEsoteric)" << std::endl;
        sourceFileLoader << "1" << std::endl;
	sourceFileLoader << "Debuggin_Level:(0-none,1-some,2-lots)" << std::endl;
	sourceFileLoader << "0" <<std::endl;
	sourceFileLoader << "Results_Output_Type:(0-none,1-to_Screen_summary,2-to_Screen_Detailed,3-to_File,4-to_File_detailed)" << std::endl;
	sourceFileLoader << "1" <<std::endl;
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
