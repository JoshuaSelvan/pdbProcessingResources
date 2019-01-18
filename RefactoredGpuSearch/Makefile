FLAGS=-g
OBJS=AtomNameToNumHashTable.o cpuBruteForceRangeSearch.o dataStructures.o ESBTL_ProteinHandler.o gpuBruteForceRangeSearch.o gpuKdTreeRangeSearch.p miscFunctions.o rangeSearch.o SearchStructureConstruction.o


gpuRangeSearchRunner: ${OBJS} kernel.cu
	nvcc ${FLAGS} ${OBJS} kernel.cu -o gpuRangeSearchRunner 
clean:
	rm -f COMPILE*
	rm -f *.o
	rm -f gpuRangeSearchRunner 

ESBTL_ProteinHandler.o: ESBTL_ProteinHandler.cpp
	g++ -c ${FLAGS} -o ESBTL_ProteinHandler.o  -I /home/jselvan/ESBTL-1.0-beta01/include/ ESBTL_ProteinHandler.cpp


%.o:	%.cpp %.h
	g++ -c ${FLAGS} -o $@ $<

%.o:	%.cu %.cuh
	nvcc -c ${FLAGS} -o $@ $<
