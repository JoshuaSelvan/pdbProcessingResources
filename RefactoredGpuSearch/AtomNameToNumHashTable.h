#ifndef ATOMHASHTABLE
#define ATOMHASHTABLE

#include <string>
#include <iostream>
#include <vector>



class AtomToNumHashTable
{
		private:
			struct HashRow
			{
				std::vector<std::string> Pos;
				//std::string Pos[9000];
			};
		 HashRow hashTable[27];
			
     
     public:
     
     AtomToNumHashTable(); //initialise class by sizing the private hashTable
     short retrieveHashValue(std::string atomName); 
    
    
    
    
};


#endif