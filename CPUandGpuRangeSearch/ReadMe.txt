This project contains code to run brute force and kd-tree range searches on unzipped PDB files with either a CPU or GPU inclusive algorithms.

If you have downloaded or recieved this project from elsewhere, you may have to download some standard libraries to get the makefile to run. This project requires the c++ library boost to be installed on the system and was built was CUDA V 10.1 with the g++ compiler V 4.8.5. This project also makes use of the ESBTL thrd party library but that should already be present within a "thirdParty" folder within the project.

With the above resources, the project should build by simply using the "make" command within the folder - which will run the compile commands in the "Makefile"

Once compiled, the program can be run by typing in ./gpuRangeSearchRunner.

The program takes its settings from a file called runDetails.txt.
This project includes a sample input list and file of input PDB files.To select your own target files you must create a list listing the locaitons of all the PDB files you wish to process and then give the path to that file in the runDetails.txt.
You may select the algorithm you wish to run from within the runDetails.txt by selecting a number in the demarkated line. The main algorithms are numbers 1,2,5,8.

In the runDetails file you will need to specify:
 - The location of your PDB location list
 - Amount of debugging
 - Results output Type.
 - The type of range search to be run
 - The atom types being sought
 - The maximum distance between the atoms for results
 - The amount of space to allocate for each length of protein being processed and how long those lengths shall be
