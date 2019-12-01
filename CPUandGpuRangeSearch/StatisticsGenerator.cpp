#include "StatisticsGenerator.h"


struct atomCoordsLocal
{
	int xyz[3];
};



int squaredDistance3DLocal(atomCoords targetAtom, atomCoords currentAtom)
{
	int temp = pow((targetAtom.xyz[0] - currentAtom.xyz[0]), 2) + pow((targetAtom.xyz[1] - currentAtom.xyz[1]), 2) + pow((targetAtom.xyz[2] - currentAtom.xyz[2]), 2);
	return temp;
}

void generateStatistics(rangeSearchSettings& settings, ProteinDataHandler heldProteinSets, AtomToNumHashTable atomReferenceTable)
{
	
	int sumOfAtomTypesAcrossAllElements = 0;
	for (int currentEntry = 0; currentEntry < heldProteinSets.ProteinDataHolder[0].heldEntries; currentEntry++)
	{
			sumOfAtomTypesAcrossAllElements = sumOfAtomTypesAcrossAllElements + heldProteinSets.ProteinDataHolder[0].compositionCountsList[currentEntry];
		
		
	
	}
	double averageElementTypeCount = sumOfAtomTypesAcrossAllElements / heldProteinSets.ProteinDataHolder[0].heldEntries;

	int currentAtomTypeCount=0;
	double sumOfVariences = 0;
	for (int currentEntry = 0; currentEntry < heldProteinSets.ProteinDataHolder[0].heldEntries; currentEntry++)
	{
				currentAtomTypeCount = heldProteinSets.ProteinDataHolder[0].compositionCountsList[currentEntry];
				sumOfVariences = sumOfVariences + (currentAtomTypeCount - averageElementTypeCount)*(currentAtomTypeCount - averageElementTypeCount);



	}
	double averageVarience = sumOfVariences /  heldProteinSets.ProteinDataHolder[0].heldEntries;






	//check the average number of 4 proximity neighbours that every atom in a protein has

	//check the varience in the number of 4 proximity neighbours that every atom in a protein has
	double sumOfAllAveragedNeigbhourCounts=0; 
	for (int currentEntry = 0; currentEntry < heldProteinSets.ProteinDataHolder[0].heldEntries; currentEntry++)
	{
		int totalNeighboursInProtein=0;
		for (int currentAtom = 0; currentAtom < heldProteinSets.ProteinDataHolder[0].proteinLengthCounts[currentEntry]; currentAtom++)
		{
			int totalNeighboursOfAtom = 0;


			atomCoords thisAtom;
			atomCoords otherAtom;

			setAtomCoords(thisAtom, heldProteinSets.ProteinDataHolder[0].xCoordsSets[currentAtom + heldProteinSets.ProteinDataHolder[0].MaxEntrySize*currentEntry], heldProteinSets.ProteinDataHolder[0].yCoordsSets[currentAtom + heldProteinSets.ProteinDataHolder[0].MaxEntrySize*currentEntry], heldProteinSets.ProteinDataHolder[0].zCoordsSets[currentAtom + heldProteinSets.ProteinDataHolder[0].MaxEntrySize*currentEntry]);



			
			for (int OtherAtom = 0; OtherAtom < heldProteinSets.ProteinDataHolder[0].proteinLengthCounts[currentEntry]; OtherAtom++)
			{
				setAtomCoords(otherAtom, heldProteinSets.ProteinDataHolder[0].xCoordsSets[OtherAtom + heldProteinSets.ProteinDataHolder[0].MaxEntrySize*currentEntry], heldProteinSets.ProteinDataHolder[0].yCoordsSets[OtherAtom + heldProteinSets.ProteinDataHolder[0].MaxEntrySize*currentEntry], heldProteinSets.ProteinDataHolder[0].zCoordsSets[OtherAtom + heldProteinSets.ProteinDataHolder[0].MaxEntrySize*currentEntry]);

				if (squaredDistance3DLocal(thisAtom, otherAtom)<=(4000*4000))
				{
					totalNeighboursOfAtom++;
				}
			}

			totalNeighboursInProtein = totalNeighboursInProtein + totalNeighboursOfAtom - 1;
		}

		sumOfAllAveragedNeigbhourCounts = sumOfAllAveragedNeigbhourCounts + totalNeighboursInProtein / heldProteinSets.ProteinDataHolder[0].proteinLengthCounts[currentEntry];

	}

	double averageNeighborCount = sumOfAllAveragedNeigbhourCounts / heldProteinSets.ProteinDataHolder[0].heldEntries;
	
	
	std::cout << "Average atom type: " << averageElementTypeCount << std::endl;
	std::cout << "Average atom type varience:" << averageVarience << std::endl;
	std::cout << "Average neigbours of atoms in proteins: " << averageNeighborCount << std::endl;

	//check the average dimensions of proteins in a set

	//check the varience in dimensions in proteins in a set

	std::cout << "fish" << std::endl;

	return;

}