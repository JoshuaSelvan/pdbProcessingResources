
#include "AtomNameToNumHashTable.h"




AtomToNumHashTable::AtomToNumHashTable()
{
	
	for (int i = 0; i < 27; i++)
	{
		hashTable[i].Pos.resize(9001);
		
		for (int j = 0; j < 9000; j++)
			AtomToNumHashTable::hashTable[i].Pos[j] = "Null";
	}
};



short AtomToNumHashTable::retrieveHashValue(std::string atomName)
{
        char firstLetterOfAtom = tolower(atomName[0]);
        short currentPos = 0;
        short storageRow = firstLetterOfAtom - 97;

        bool atomAlreadyStored = true;
        if ((storageRow>-1) && (storageRow<26))
        {
                while (atomAlreadyStored)
                {
					if (AtomToNumHashTable::hashTable[storageRow].Pos[currentPos] == atomName)
                        {
                                return currentPos * 100 + storageRow;
                        }
                        else if (hashTable[storageRow].Pos[currentPos] == "Null")
                        {
                                //std::cout<<"Atom was not previously in hashTable"<<std::endl; uncomment to see when atom names are being added to the hashtable
                                hashTable[storageRow].Pos[currentPos] = atomName;
                                return currentPos * 100 + storageRow;
                        }
                        else
                                currentPos++;
                }
        }
        else
        {
                while (atomAlreadyStored)
                {
                        if (hashTable[26].Pos[currentPos] == atomName)
                        {
                                return currentPos * 100 + storageRow;
                        }
                        else if (hashTable[26].Pos[currentPos] == "Null")
                        {
                                hashTable[26].Pos[currentPos] = atomName;
                                return currentPos * 100 + storageRow;
                        }
                        else
                                currentPos++;
                }
        }
        return 0;
}