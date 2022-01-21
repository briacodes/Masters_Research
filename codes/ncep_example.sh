#!/bin/sh
#SBATCH -v
#SBATCH --job-name=ncep
#SBATCH --mem-per-cpu=7000
#SBATCH --partition=all
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=0:30:00
#SBATCH --output=/home/snebuda/cris/ncep_cloudtest/jobout/ncep.%j.out

exe=/apollo/jung/bandersen/ncep_cloudtest/ncep.exe
file_satinfo=/apollo/jung/bandersen/ncep_cloudtest/global_satinfo.txt

source /etc/bashrc
module purge
module load license_intel
module load intel
module load hdf
module load hdf5
module load netcdf4

input_path=/apollo/jung/snebuda/brianne/2021
output_path=/apollo/jung/bandersen/ncep_data_2021_07_30

cd ${input_path}

#file_bt=cris.431.n20.20210410_000000.20210410_003000.bt.nc
file_bt=cris.431.n20.20210730_000000.20210730_003000.bt.nc # template substitution here

froot=${file_bt:: -5}
#echo ${froot}

file_state=${froot}state.nc
file_te=${froot}te.nc
file_sens=${froot}sens.nc
#file_out=${froot}ncep.nobias.nc

file_out=${output_path}/${froot}ncep.nc

echo file out is ${file_out}

srun --cpu_bind=core --distribution=block:block ${exe} \
${file_satinfo} \
${file_bt} \
${file_state} \
${file_te} \
${file_sens} \
${file_out}

sacct -j $SLURM_JOB_ID --format=JobID,JobName,Elapsed,CPUTime,MaxRSS,NodeList%30

exit
~
