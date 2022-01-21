#!/bin/sh
#SBATCH --job-name=ncep_dr_cd_fov
#SBATCH --partition=all
#SBATCH --mem-per-cpu=60000
#SBATCH --ntasks=1
#SBATCH --time=05:00:00
#SBATCH --output=/home/bandersen/iris-home/%j.out
#
#DONT FORGET
#the directory for the batch output (above) has to exist or the job won't run
cd /home/bandersen/iris-home
module load license_intel intel hdf hdf5 netcdf4 miniconda

python all_algorithms/2020_10_10_co2/09_10.py

#srun executable if it was set up for parallel processing
#otherwise just running the python like at the command line
#will have the resources of the batch job you requested
#at the top of the script


