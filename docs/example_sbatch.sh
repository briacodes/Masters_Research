#!/bin/sh
#SBATCH --job-name=myjob
#SBATCH --partition=mars
#SBATCH --ntasks=1 
#SBATCH --time=00:10:00
#SBATCH --output=/scratch/bandersen/data/myjob.%j.out
#
cd /scratch/bandersen/data

#the directory for the batch output has to exist or the job won't run
#/scratch/bandersen/data/ is listed above in the output argument to SBATCH

echo Hello!

module load license_intel intel hdf hdf5 netcdf4 miniconda

#srun executable






