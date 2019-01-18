
#include "ESBTL_ProteinHandler.h"
#include <ESBTL/default.h>
#include "miscFunctions.h"



typedef ESBTL::Accept_none_occupancy_policy<ESBTL::PDB::Line_format<> > Accept_none_occupancy_policy;



ProteinDataHandler::ProteinDataHandler(rangeSearchSettings currentSettings)
{
    for(int i=0;i<5;i++)
    {
        ProteinDataHolder[i].heldEntries=0;
	ProteinDataHolder[i].MaxEntrySize = currentSettings.maxProteinLengthPerRange[i];
        ProteinDataHolder[i].MaxNumOfEntries=currentSettings.maxNumOfEntriesPerRange[i];
	ProteinDataHolder[i].KdTreeSize = currentSettings.maxProteinLengthPerRange[i] * 2;
        ProteinDataHolder[i].namesSets=(short*)malloc(ProteinDataHolder[i].MaxEntrySize * ProteinDataHolder[i].MaxNumOfEntries * sizeof(short));
        ProteinDataHolder[i].xCoordsSets= (int*)malloc(ProteinDataHolder[i].MaxEntrySize * ProteinDataHolder[i].MaxNumOfEntries * sizeof(int));
        ProteinDataHolder[i].yCoordsSets= (int*)malloc(ProteinDataHolder[i].MaxEntrySize * ProteinDataHolder[i].MaxNumOfEntries * sizeof(int));
        ProteinDataHolder[i].zCoordsSets= (int*)malloc(ProteinDataHolder[i].MaxEntrySize * ProteinDataHolder[i].MaxNumOfEntries * sizeof(int));
        ProteinDataHolder[i].proteinLengthCounts=(int*)malloc(ProteinDataHolder[i].MaxNumOfEntries*sizeof(int));
		ProteinDataHolder[i].kdTrees = (int*)malloc(ProteinDataHolder[i].KdTreeSize * ProteinDataHolder[i].MaxNumOfEntries * sizeof(int));
    }
    
    for(int i=0;i<5;i++)
    {
        ProteinDataHolder[i].compositionLists= (int*)malloc(ProteinDataHolder[i].MaxEntrySize * ProteinDataHolder[i].MaxNumOfEntries * sizeof(int)); 
		ProteinDataHolder[i].compositionPointers = (int*)malloc(ProteinDataHolder[i].MaxEntrySize * ProteinDataHolder[i].MaxNumOfEntries * sizeof(int));
		ProteinDataHolder[i].compositionCountsList = (int*)malloc(ProteinDataHolder[i].MaxNumOfEntries * sizeof(int));
       
    }



    
    for (int i = 0; i < 5; i++)
	{
		assert(ProteinDataHolder[i].namesSets != NULL);
		assert(ProteinDataHolder[i].xCoordsSets != NULL);
		assert(ProteinDataHolder[i].yCoordsSets != NULL);
		assert(ProteinDataHolder[i].zCoordsSets != NULL);
		assert(ProteinDataHolder[i].namesSets != NULL);
		assert(ProteinDataHolder[i].kdTrees != NULL);

		assert(ProteinDataHolder[i].compositionLists != NULL);
		assert(ProteinDataHolder[i].compositionCountsList != NULL);
		assert(ProteinDataHolder[i].proteinLengthCounts != NULL);
	}
    
	int x = 5;
    for(int i=0;i<5;i++)
    {
        setShortArrayValues(ProteinDataHolder[i].namesSets, -1, ProteinDataHolder[i].MaxEntrySize * ProteinDataHolder[i].MaxNumOfEntries);
		setIntArrayValues(ProteinDataHolder[i].xCoordsSets, 999999, ProteinDataHolder[i].MaxEntrySize * ProteinDataHolder[i].MaxNumOfEntries);
		setIntArrayValues(ProteinDataHolder[i].yCoordsSets, 999999, ProteinDataHolder[i].MaxEntrySize * ProteinDataHolder[i].MaxNumOfEntries);
		setIntArrayValues(ProteinDataHolder[i].zCoordsSets, 999999, ProteinDataHolder[i].MaxEntrySize * ProteinDataHolder[i].MaxNumOfEntries);
		setIntArrayValues( ProteinDataHolder[i].kdTrees, -1, ProteinDataHolder[i].KdTreeSize*ProteinDataHolder[i].MaxNumOfEntries); 
    }
    

}

void ProteinDataHandler::loadAllProteinsToArrays(std::string inputFileLocationsList, AtomToNumHashTable &atomReferenceTable)
{
    std::ifstream fileLoader;
    fileLoader.open(inputFileLocationsList.c_str());
	std::string currentFile;
    
    short filesLoaded =0;
    //int filesLoadedPerSet[5];
    for (int i=0;i<5;i++)
       
    ProteinDataHolder[i].heldEntries;
    
    while ((fileLoader >> currentFile)&&(filesLoaded<1000))
	{
		filesLoaded = filesLoaded + loadSingleProteinDetailsIntoArrays(currentFile, atomReferenceTable);
    }
    std::cout<<"number of files succesfully loaded into the holding arrays: "<<filesLoaded<<std::endl<<std::endl;
    for (int i=0; i<5;i++)
    {
		std::cout << "Number of entries in the " << ProteinDataHolder[i].MaxEntrySize << " range is: " << ProteinDataHolder[i].heldEntries << std::endl;
    }
    
    
};

int ProteinDataHandler::loadSingleProteinDetailsIntoArrays(std::string proteinFileLocation, AtomToNumHashTable &atomReferenceTable)
{
    int numOfAtomsInCurrentProtein = 0;
    int destinationSet = 0;
    int arrayInsertionStartPoint =0;
    
	std::cout << "Processing file: " << proteinFileLocation << std::endl;
    ESBTL::PDB_line_selector_two_systems sel;
    std::vector<ESBTL::Default_system> systems;
    ESBTL::All_atom_system_builder<ESBTL::Default_system> builder(systems, sel.max_nb_systems());
	if (ESBTL::read_a_pdb_file(proteinFileLocation, sel, builder, Accept_none_occupancy_policy()))
    {
           if (systems.empty() || systems[0].has_no_model()){
			   std::cerr << "No atoms found in file " << proteinFileLocation << std::endl;
                    return 0;
            }
            else
            {
                    for (ESBTL::Default_system::Models_iterator it_model = systems[0].models_begin(); it_model != systems[0].models_end(); ++it_model)
                    {
                           const ESBTL::Default_system::Model& model = *it_model;
                           for (ESBTL::Default_system::Model::Atoms_const_iterator it_atm = model.atoms_begin(); it_atm != model.atoms_end(); ++it_atm)
                           {
                                   numOfAtomsInCurrentProtein++;
                           }
                    } 
                    
                    
                    if(numOfAtomsInCurrentProtein!=0)
                    {
						int t = 0;
                        for ( t = 0; t < 5; t++)
			            {
				           if (numOfAtomsInCurrentProtein < ProteinDataHolder[t].MaxEntrySize)
				           {
					           destinationSet = t;
					           t = t + 100;
				           }
			           }
                       if((t>99)&&(ProteinDataHolder[destinationSet].heldEntries<ProteinDataHolder[destinationSet].MaxNumOfEntries))
                       { 
                           int localCounter=0;
                           arrayInsertionStartPoint= ProteinDataHolder[destinationSet].heldEntries* ProteinDataHolder[destinationSet].MaxEntrySize;
                           for (ESBTL::Default_system::Models_iterator it_model = systems[0].models_begin(); it_model != systems[0].models_end(); ++it_model)
                           {
                               const ESBTL::Default_system::Model& model = *it_model;
                               for (ESBTL::Default_system::Model::Atoms_const_iterator it_atm = model.atoms_begin(); it_atm != model.atoms_end(); ++it_atm)
                               { 
								   ProteinDataHolder[destinationSet].namesSets[arrayInsertionStartPoint + localCounter] = atomReferenceTable.retrieveHashValue(it_atm->atom_name());
                                    ProteinDataHolder[destinationSet].xCoordsSets[arrayInsertionStartPoint+localCounter] = (it_atm->x()*1000);
                                    ProteinDataHolder[destinationSet].yCoordsSets[arrayInsertionStartPoint+localCounter] = (it_atm->y()*1000);
                                    ProteinDataHolder[destinationSet].zCoordsSets[arrayInsertionStartPoint+localCounter] = (it_atm->z()*1000);
                                    localCounter++;
                               }
                           }
                           ProteinDataHolder[destinationSet].proteinLengthCounts[ProteinDataHolder[destinationSet].heldEntries]=localCounter;
                           ProteinDataHolder[destinationSet].heldEntries++;
                           return 1;
                       }
                       else
                       {
                           std::cout<<"file did not fall in valid range"<<std::endl;
                           return 0; 
                       }
                    }
            }
    }
    else
    {
        std::cout<<"error occured with loading file"<<std::endl;
        return 0;
    }


	return -3;
    
};

int checkNumberOfProteinFilesInList(std::string inputFileLocationList)
{
	int numberOfFiles=0;

	std::ifstream fileLoader;
	fileLoader.open(inputFileLocationList.c_str());
	std::string currentFile;

	while (fileLoader >> currentFile)
	{
		numberOfFiles++;
	}
	std::cout << "Total number of files is current set" << numberOfFiles << std::endl;
	return numberOfFiles;
};







