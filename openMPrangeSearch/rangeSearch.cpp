//#include <mpi.h>
#include <omp.h>
#include <stdlib.h>
#include <iostream>
#include <string>

#include <stdio.h>
#include <time.h>
#include "ESBTL_ProteinHandler.h"

#include "miscFunctions.h"
#include "dataStructures.h"
#include <string>

#include "cpuBruteForceSearch.h"
#include "cpuKdTreeSearch.h"
#include <chrono>
using namespace std;
using namespace std::chrono;



int main(int argc, char** argv) {
//Main program
 typedef std::chrono::high_resolution_clock Time;
 typedef std::chrono::milliseconds ms;
 typedef std::chrono::duration<float> fsec;


 clock_t start, stop;

 high_resolution_clock::time_point t1;
 high_resolution_clock::time_point t2; 
 rangeSearchSettings currentSettings;
 std::string runSettingsFile = "runDetails.txt";
 loadRunDetails(runSettingsFile, currentSettings);
 start = clock();


std::cout<<"Number of files is: "<<currentSettings.numberOfFiles<<std::endl;
AtomToNumHashTable atomReferenceTable;
ProteinDataHandler heldProteinSets(currentSettings);


std::ifstream fileLoader;
fileLoader.open(currentSettings.inputListFileLocation.c_str());
vector<std::string> inputLocationsList;

inputLocationsList.resize(currentSettings.numberOfFiles);
for(int i=0;i<currentSettings.numberOfFiles;i++)
	fileLoader>>inputLocationsList[i];

t1 = high_resolution_clock::now();

//Parallel file loading approach
/*#pragma omp parallel for shared(heldProteinSets) 
for(int i=0;i<currentSettings.numberOfFiles;i++)
{
//Load 1 file into memory
std::string fileToBeProcessed =inputLocationsList[i]; 
std::cout<<"about to process entry "<<i<<std::endl;
heldProteinSets.loadSingleProteinToArraysOmp(atomReferenceTable,currentSettings , i, fileToBeProcessed);
}
*/

//Linear file loading approach
heldProteinSets.loadAllProteinsToArrays(currentSettings.inputListFileLocation, atomReferenceTable, currentSettings);
t2 = high_resolution_clock::now();
print_chrono_elapsed(t1,t2,"Load from esbtl CHRONOVERSION: ");

t1 = high_resolution_clock::now();
//OMP_Code_loading
{
heldProteinSets.nullSecondarySearchStructures();
heldProteinSets.FormatSecondaryPositionStructuresOMP();
}

//non_OMP_codes
//heldProteinSets.setSecondarySearchStructure();
t2 = high_resolution_clock::now();
print_chrono_elapsed(t1,t2,"CPU SECONDARY SEARCH STRUCTURE CONSTRUCTION TIME: ");



t1 = high_resolution_clock::now();
constructKdTreesOnLoadedDataOnCPUOMP(heldProteinSets);
t2 = high_resolution_clock::now();
 print_chrono_elapsed(t1,t2,"kd Tree Construction chronotime: ");

//cpuKdTreeRangeSearch(currentSettings, heldProteinSets, atomReferenceTable);
//heldProteinSets.DisplayHeldEntriesPerRange();
//start=clock();

//This runs a normal sequential run
//cpuBruteForceRangeSearchAllLoadedSets(currentSettings, heldProteinSets, atomReferenceTable);
//return 0;

//This loops is mearly to repeat the experiment 5 times.
//
//

//cout<< "PERFORMING SECONDARY STRUCTURE SEARCH"<<endl;
//for(int r=0; r<0; r++)
//{

/*
if(r==0)
{
currentSettings.AtomTypeOne="CB";
currentSettings.AtomTypeTwo="CG";
}
if(r==1)
{
currentSettings.AtomTypeOne="C1'";
currentSettings.AtomTypeTwo="C2'";
}
if(r==2)
{
currentSettings.AtomTypeOne="rrr";
currentSettings.AtomTypeTwo="ttt";
}
*/
if (currentSettings.searchType == 1)
{
cout<< "PERFORMING SECONDARY STRUCTURE SEARCH"<<endl;
cout<<"Processing Atom Pair: "<<currentSettings.AtomTypeOne<<" + "<<currentSettings.AtomTypeTwo<<endl;
t1 = high_resolution_clock::now();
for(int j=0;j<5;j++)
{
        std::cout << std::endl;
        std::cout << "Processing Range set: " << j << std::endl;
        std::cout << "Number of present entries is: " << heldProteinSets.ProteinDataHolder[j].heldEntries << std::endl;
#pragma omp parallel for 
for(int i=0;i<heldProteinSets.ProteinDataHolder[j].heldEntries; i++)
{
processSingleFile(currentSettings, heldProteinSets, atomReferenceTable, j, i);
}

}

t2 = high_resolution_clock::now();
print_chrono_elapsed(t1,t2,"elapsed chronotime: ");
std::cout << std::endl;
}
//}

std::cout << std::endl;
std::cout << std::endl;
std::cout << std::endl;

//return 0;

//for (int r=2; r<3; r++)
//{
if (currentSettings.searchType ==2)
{
t1 = high_resolution_clock::now();

 std::cout << "PERFORMING KD TREE RANGE SEARCH" << std::endl;
//if(r==0)
//{
//currentSettings.AtomTypeOne="CB";
//currentSettings.AtomTypeTwo="CG";
//}
//if(r==1)
//{
//currentSettings.AtomTypeOne="OG";
//currentSettings.AtomTypeTwo="N";
//}
//if(r==2)
//{
///currentSettings.AtomTypeOne="rrr";
//currentSettings.AtomTypeTwo="ttt";
//}


cout<<"Processing Atom Pair: "<<currentSettings.AtomTypeOne<<" + "<<currentSettings.AtomTypeTwo<<endl;

for (int i = 0; i < 5; i++)
{
	if (heldProteinSets.ProteinDataHolder[i].heldEntries>0)
	{
		std::cout << std::endl;
		std::cout << "Processing Range set: " << i << std::endl;
		std::cout << "Number of present entries is: " << heldProteinSets.ProteinDataHolder[i].heldEntries << std::endl;
        	#pragma omp parallel for
		for (int currentEntry = 0; currentEntry < heldProteinSets.ProteinDataHolder[i].heldEntries; currentEntry++)
		{
	                processSingleKdTreeRangeSearch(currentSettings,  heldProteinSets,  atomReferenceTable, i, currentEntry);

        	}
        }
}


t2 = high_resolution_clock::now();


 print_chrono_elapsed(t1,t2,"elapsed kdTree chronoTime: ");
//	print_elapsed(start2, end2, "BETTER TIMER = = = = = Total run time for kdTree: ");
//	delta = ((end2.tv_sec  - start2.tv_sec) * 1000000u + end2.tv_usec - start2.tv_usec) / 1.e6;
//        cout<<"Better time = "<< delta<<endl;
std::cout << std::endl;
}

//}
	return 0;


//Example2
/*	int currentNum=0;
	int* list = (int*)malloc(sizeof(int)*12);
	for (int i=0;i<12;i++)
        	list[i]=9999; 
     #pragma omp parallel shared(list) num_threads(6)
        {
        int ID = omp_get_thread_num();
	//char numstr [10];
        //string outputFileName = "InputListForID";
        //sprintf(numstr, "%d", ID);
        //outputFileName += numstr;
        //outputFileName += ".txt";
	//list[ID]=ID;
	//ofstream myfile;
  	//myfile.open (outputFileName);
  	//myfile <<"bla \n";
  	int localVal;
  	#pragma omp atomic capture
	localVal= currentNum++;
  	//myfile <<list[ID];
  	list[ID]=localVal;
  	//myfile<<localVal;
	//myfile.close();
	}
ofstream lastFile;
lastFile.open("fullOutput.txt");
 for (int i=0;i<12;i++)
               lastFile<<list[i]<<" \n";
lastFile.close();
return 0;

*/

//Example1
/*int nn =4;
int* list = (int*)malloc(sizeof(int)*12); 
for (int i=0;i<12;i++)
	list[i]=i;
for (int i=0;i<12;i++)
        cout<<list[i]<<endl;
        


	#pragma omp parallel shared(list) num_threads(6)
	{
	int ID = omp_get_thread_num();
	int Total=0;
	//cout<<"ID is set to: "<<ID<<endl;
	Total=list[ID]+list[ID+6];
	cout<<"Thread "<<ID<<" has been initialised "<<std::endl;
	list[ID] =Total;	
	


	}
for(int i=0;i<12;i++)
	cout<<list[i]<<endl;
*/





//        int totalNumberOfFiles = checkNumberOfProteinFilesInList(RangeSearchSettings.inputListFileLocation);
//	std::cout<<totalNumberOfFiles<<endl;

	return 0;
}
