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

t2 = high_resolution_clock::now();
print_chrono_elapsed(t1,t2,"CPU SECONDARY SEARCH STRUCTURE CONSTRUCTION TIME: ");



t1 = high_resolution_clock::now();
constructKdTreesOnLoadedDataOnCPUOMP(heldProteinSets);
t2 = high_resolution_clock::now();
 print_chrono_elapsed(t1,t2,"kd Tree Construction chronotime: ");



cout<< "PERFORMING SECONDARY STRUCTURE SEARCH"<<endl;
for(int r=0; r<0; r++)
{

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


cout<<"Processing Atom Pair: "<<currentSettings.AtomTypeOne<<" + "<<currentSettings.AtomTypeTwo<<endl;
t1 = high_resolution_clock::now();
for(int j=0;j<5;j++)
{
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

std::cout << std::endl;
std::cout << std::endl;
std::cout << std::endl;

for (int r=2; r<3; r++)
{
t1 = high_resolution_clock::now();

 std::cout << "PERFORMING KD TREE RANGE SEARCH" << std::endl;
if(r==0)
{
currentSettings.AtomTypeOne="CB";
currentSettings.AtomTypeTwo="CG";
}
if(r==1)
{
currentSettings.AtomTypeOne="OG";
currentSettings.AtomTypeTwo="N";
}
if(r==2)
{
currentSettings.AtomTypeOne="rrr";
currentSettings.AtomTypeTwo="ttt";
}


cout<<"Processing Atom Pair: "<<currentSettings.AtomTypeOne<<" + "<<currentSettings.AtomTypeTwo<<endl;

for (int i = 0; i < 5; i++)
{
	if (heldProteinSets.ProteinDataHolder[i].heldEntries>0)
	{
        	#pragma omp parallel for
		for (int currentEntry = 0; currentEntry < heldProteinSets.ProteinDataHolder[i].heldEntries; currentEntry++)
		{
	                processSingleKdTreeRangeSearch(currentSettings,  heldProteinSets,  atomReferenceTable, i, currentEntry);

        	}
        }
}


t2 = high_resolution_clock::now();


 print_chrono_elapsed(t1,t2,"elapsed kdTree chronoTime: ");
std::cout << std::endl;


}
	return 0;

	return 0;
}
