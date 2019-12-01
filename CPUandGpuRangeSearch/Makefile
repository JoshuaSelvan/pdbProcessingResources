FLAGS=-g 
OBJS=AtomNameToNumHashTable.o cpuBruteForceSearch.o cpuKdTreeSearch.o dataStructures.o ESBTL_ProteinHandler.o gpuBruteForceSearch.o gpuKdTreeRangeSearch.o interface.o miscFunctions.o rangeSearch.o SearchStructureConstruction.o printToFileHandler.o gpuMemoryLimitTester.o StatisticsGenerator.o constructKdTreesOnLoadedDataWithGPU.o


guRangeSearchRunner: ${OBJS} kernel.cu
	nvcc ${FLAGS} ${OBJS} kernel.cu -std=c++11 -o gpuRangeSearchRunner  
clean:
	rm -f COMPILE*
	rm -f *.o
	rm -f gpuRangeSearchRunner

ESBTL_ProteinHandler.o: ESBTL_ProteinHandler.cpp
	g++ -static-libstdc++ -c ${FLAGS} -o ESBTL_ProteinHandler.o -I /home/jselvan/ESBTL-1.0-beta01/include/ ESBTL_ProteinHandler.cpp -std=c++11   


%.o:    %.cpp %.h
	/usr/local/cuda-10.1/bin/nvcc -std=c++11 -c ${FLAGS} -o $@ $<

%.o:    %.cu %.cuh
	nvcc -std=c++11 -c ${FLAGS} -o $@ $<
