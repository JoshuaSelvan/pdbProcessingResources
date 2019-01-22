#ifndef DATASTRUCTURES
#define DATASTRUCTURES
#include <string>
#include <vector>





class rangeSearchSettings
{
public:
	int maxNumOfEntriesPerRange[5];
	int maxProteinLengthPerRange[5];
    std::string inputListFileLocation; 
    std::string AtomTypeOne;
    std::string AtomTypeTwo;
    int requiredProximity;
    int numberOfFiles;
	int searchType;
	int debugLevel;
	int resultsPrintFormat;

	int multiBatchRangeSizes[5];
};

struct multipleProteinSet
    {
        int heldEntries; //The number of entries currently being stored in the arrays
	    int MaxEntrySize; //The maximum length of protiens stored in this set.
        int KdTreeSize; //The set size of the kdTrees associated with these lengths of proteins
	    int MaxNumOfEntries; //Max number of records this array is set to hold
    
        //The attributes of each atom in each protein are stored in these - with consecutive proteins starting at set intervals
        int *xCoordsSets;
	    int *yCoordsSets;
	    int *zCoordsSets;
	    short *namesSets;
        int *kdTrees; //The kdTrees mapping to the above storage arrays
	    int*proteinLengthCounts; //The length of every protein chain stored in the above arrays
        
        
        //These are only used if the range search utilises a secondary structure for locating all instances of a set of atoms instantly
		//Secondary reference arrays.
	    int*compositionCountsList; //Contains the number of unique atoms held in each protein chain
	    int*compositionLists; // contains the name of each atom type in the protein followed by how many instances of the atom are present.
	    int*compositionPointers; //Contains pointers to where all secondary elements are from in the protein (but reorganises them so that the same atom types are placed together in storage
    
		int* xyzCoordsSets;

    }; 
    
    
    
    struct rangeSearchArrays
    {
	int* atomAPositionList;
	int* atomACount;

	int* MatchesCount;
	int* nextSearchCount;
	int* completionFlag;

	int* atomACurrentSearchDimensions;
	int* atomACurrentSearchKdTreePositions;

	int* atomAMatches;
	int* atomBMatches;

    
	int blocks;
	int threads;
	//Test arrays which play no active role in the code
	int* atomBPositionList;

    };

	struct gpuBruteForceSingleEntryResources
	{
		int* h_atomAPositionList;
		int* h_atomACount;

		int* h_atomAMatches;
		int* h_atomBMatches;
		int* h_MatchesCount;
		int* h_nextSearchCount;
		int* h_completionFlag;

		int* d_atomAPositionList;
		int* d_atomACurrentSearchDimensions;
		int* d_atomACurrentSearchKdTreePositions;
		int* d_atomACount;

		int* d_atomAMatches;
		int* d_atomBMatches;
		int* d_MatchesCount;

		int* d_nextSearchCount;
		int* d_completionFlag;

		int blocks;
		int threads;

		//used for kd-tree search
		int *d_xyzCoordsSets;

		//used for brute=force search
		int *d_xCoords;
		int *d_yCoords;
		int *d_zCoords;


		short *d_names;
		int *d_kdTree;
		////////////////////////////////////////////////////////////
		//processing resources


		//Test arrays which play no active role in the code
		//int* h_atomBPositionList;
		//int* d_atomBPositionList;
		//int* h_atomACurrentSearchDimensions;
		//int* h_atomACurrentSearchKdTreePositions;

		int concurrentThreads;

		int* h_resultsCount;
		int* h_resultsA;
		int* h_resultsB;

		int *h_aCount;
		int *h_bCount;


		int* h_elementAList;  //I used these for loading the element arrays back from the device to check what was in them, otherwise commented out.
		int* h_elementBList;

		int* d_elementAList; //I used these for loading the element arrays back from the device to check what was in them, otherwise commented out.
		int* d_elementBList;

		int* d_resultsCount;
		int* d_resultsA;
		int* d_resultsB;

		int* d_aCount;
		int* d_bCount;
	};

	struct gpuRangeSearchResources{
		int* h_atomAPositionList;
		int* h_atomACount;

		int* h_atomAMatches;
		int* h_atomBMatches;
		int* h_MatchesCount;
		int* h_nextSearchCount;
		int* h_completionFlag;

		int* d_atomAPositionList;
		int* d_atomACurrentSearchDimensions;
		int* d_atomACurrentSearchKdTreePositions;
		int* d_atomACount;

		int* d_atomAMatches;
		int* d_atomBMatches;
		int* d_MatchesCount;

		int* d_nextSearchCount;
		int* d_completionFlag;

		int blocks;
		int threads;

		//used for kd-tree search
		int *d_xyzCoordsSets;

		//used for brute=force search
		int *d_xCoordsSets;
		int *d_yCoordsSets;
		int *d_zCoordsSets;


		short *d_namesSets;
		int *d_kdTreeSets;



		//Test arrays which play no active role in the code
		int* h_atomBPositionList;
		int* d_atomBPositionList;
		int* h_atomACurrentSearchDimensions;
		int* h_atomACurrentSearchKdTreePositions;


	};


	struct gpuBruteForceResources{
		int blocks;
		int threads;
		int concurrentThreads;

		int* h_resultsCount; 
		int* h_resultsA; 
		int* h_resultsB; 

		int *h_aCount; 
		int *h_bCount;


		int* h_elementAList;  //I used these for loading the element arrays back from the device to check what was in them, otherwise commented out.
		int* h_elementBList; 

		int* d_elementAList; //I used these for loading the element arrays back from the device to check what was in them, otherwise commented out.
		int* d_elementBList;

		int* d_resultsCount;
		int* d_resultsA;
		int* d_resultsB;

		int* d_aCount;
		int* d_bCount;
	};


struct dataLoader
{
	int *numOfAtoms;

	//used in bitonic sort
	int *xValues;
	int *yValues;
	int *zValues;

	int *kdTreeholder;

	int *xBackPositionReferenceList;
	int *yBackPositionReferenceList;
	int *zBackPositionReferenceList;

	int *xForwardPositionReferenceList;
	int *yForwardPositionReferenceList;
	int *zForwardPositionReferenceList;


};


struct dataLoaderWithGpuExtensions
{
	int *numOfAtoms;

	
	int *xValues;
	int *yValues;
	int *zValues;

	int *kdTreeholder;

	int *xBackPositionReferenceList;
	int *yBackPositionReferenceList;
	int *zBackPositionReferenceList;

	int *xForwardPositionReferenceList;
	int *yForwardPositionReferenceList;
	int *zForwardPositionReferenceList;

	int *d_xBackPositionReferenceList;
	int *d_yBackPositionReferenceList;
	int *d_zBackPositionReferenceList;

	int *d_xForwardPositionReferenceList;
	int *d_yForwardPositionReferenceList;
	int *d_zForwardPositionReferenceList;

	int *d_xCoordsHolder;
	int *d_yCoordsHolder;
	int *d_zCoordsHolder;

	int*d_kdTree;


};


struct kdTask {
	int kdDestination;
	int minX;
	int maxX;
	int minY;
	int maxY;
	int minZ;
	int maxZ;
};

struct multiRunDetailsSet{

	int proximitiesCount;
	std::vector<float> proximities;
	std::vector<int> entriesToSearchEachRound;
	int currentRange ;
	std::string atom1;
	std::string atom2;
	int searchType;
	int batchSizesPerRange[5];
};


struct atomCoords
{
    int xyz[3];   
};



void setAtomCoords(atomCoords  &coordsHolder, int X, int Y, int Z);



#endif
