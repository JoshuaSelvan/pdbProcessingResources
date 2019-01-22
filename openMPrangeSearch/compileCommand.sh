
#example command to compile code
g++ -o ompV7 -fopenmp rangeSearch.cpp miscFunctions.cpp AtomNameToNumHashTable.cpp cpuBruteForceSearch.cpp printToFileHandler.cpp dataStructures.cpp cpuKdTreeSearch.cpp SearchStructureConstruction.cpp  ESBTL_ProteinHandler.cpp -I /home/jselvan/ESBTL-1.0-beta01/include/ -std=gnu++11

#command to set the number of threads that openMP will use
export OMP_NUM_THREADS=2
export OMP_NUM_THREADS=4

