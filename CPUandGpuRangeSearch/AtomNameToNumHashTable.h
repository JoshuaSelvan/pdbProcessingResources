#ifndef ATOMHASHTABLE
#define ATOMHASHTABLE

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>   


class AtomToNumHashTable
{
		private:
			struct HashRow
			{
				std::vector<std::string> Pos;
			};
		 HashRow hashTable[27];
			
     
     public:  
     AtomToNumHashTable(); //initialise class by sizing the private hashTable
     short retrieveHashValue(std::string atomName); 
 
};


#endif