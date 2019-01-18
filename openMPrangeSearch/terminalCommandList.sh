g++ -o ompiV2 -fopenmp rangeSearch.cpp miscFunctions.cpp AtomNameToNumHashTable.cpp  dataStructures.cpp SearchStructureConstruction.cpp  ESBTL_ProteinHandler.cpp -I /home/jselvan/ESBTL-1.0-beta01/include/ -std=gnu++11







g++ -o ompiV2 -fopenmp rangeSearch.cpp miscFunctions.cpp AtomNameToNumHashTable.cpp cpuBruteForceSearch.cpp printToFileHandler.cpp dataStructures.cpp SearchStructureConstruction.cpp  ESBTL_ProteinHandler.cpp -I /home/jselvan/ESBTL-1.0-beta01/include/ -std=gnu++11



g++ -o ompiV2 -fopenmp rangeSearch.cpp miscFunctions.cpp AtomNameToNumHashTable.cpp cpuBruteForceSearch.cpp printToFileHandler.cpp dataStructures.cpp cpuKdTreeSearch.cpp SearchStructureConstruction.cpp  ESBTL_ProteinHandler.cpp -I /home/jselvan/ESBTL-1.0-beta01/include/ -std=gnu++11








for valgrind:

g++ -g -o ompiV2 -fopenmp rangeSearch.cpp miscFunctions.cpp AtomNameToNumHashTable.cpp cpuBruteForceSearch.cpp printToFileHandler.cpp dataStructures.cpp cpuKdTreeSearch.cpp SearchStructureConstruction.cpp  ESBTL_ProteinHandler.cpp -I /home/jselvan/ESBTL-1.0-beta01/include/ -std=gnu++11

