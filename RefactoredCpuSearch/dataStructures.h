#ifndef DATASTRUCTURES
#define DATASTRUCTURES
#include <string>





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

struct kdTask {
	int kdDestination;
	int minX;
	int maxX;
	int minY;
	int maxY;
	int minZ;
	int maxZ;
};




struct atomCoords
{
    int xyz[3];   
};



void setAtomCoords(atomCoords  &coordsHolder, int X, int Y, int Z);



#endif