#PBS -N ompiTest
#PBS -q medium
#PBS -l walltime=00:20:00,mem=16GB
#PBS -l nodes=1

module load openmpi-x86_64

cd /home/jselvan/RefactoredCode/openMPrangeSearch

#mpiCC rangeSearch.cpp  #mpiCC is used to compile c++ while mpicc is used to compile c
rm ompiTest.* OmpiRangeSearch3
mpiCC -std=c++11 -o OmpiRangeSearch3 rangeSearch.cpp miscFunctions.cpp cpuKdTreeSearch.cpp AtomNameToNumHashTable.cpp cpuBruteForceSearch.cpp dataStructures.cpp rangeSearchHandler.cpp SearchStructureConstruction.cpp  ESBTL_ProteinHandler.cpp -I /home/jselvan/ESBTL-1.0-beta01/include/

