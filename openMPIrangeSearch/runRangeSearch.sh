#PBS -N RUNOMPIRESULTS
#PBS -q WitsLong
#PBS -l walltime=01:30:00,mem=16GB
#PBS -l nodes=4

module load openmpi-x86_64

cd /home/jselvan/openMPDTest/rangeSearch

#echo "Check which machines we've been given"
#cat ${PBS_NODEFILE}

#mpi-start -t openmpi -np 4 a.out
mpi-start -t openmpi -np 2 --hostfile hostfiles.txt OmpiRangeSearch2 
# or mpirun -np 4 ./a.out -- which gives no errors
mpi-start -t openmpi -np 4 OmpiRangeSearch


 mpirun -np 2 --hostfile hostfiles.txt OmpiRangeSearch2

