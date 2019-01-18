FLAGS=-g
OBJS= rangeSearch.o miscFunctions.o dataStructures.o cpuKdTreeSearch.o cpuBruteForceSearch.o AtomNameToNumHashTable.o ESBTL_ProteinHandler.o SearchStructureConstruction.o


cpuRangeSearch: ${OBJS} RefactoredRangeSearchV1.cpp
	g++ ${FLAGS} ${OBJS} RefactoredRangeSearchV1.cpp -o cpuRangeSearch
clean:
	rm -f COMPILE*
	rm -f *.o
	rm -f cpuRangeSearch

ESBTL_ProteinHandler.o: ESBTL_ProteinHandler.cpp
	g++ -c ${FLAGS} -o ESBTL_ProteinHandler.o  -I /home/jselvan/ESBTL-1.0-beta01/include/ ESBTL_ProteinHandler.cpp


%.o:	%.cpp %.h
	g++ -c ${FLAGS} -o $@ $<

