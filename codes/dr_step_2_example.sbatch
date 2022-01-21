#!/bin/sh
#SBATCH --job-name=ncep__dr_cd_fov
#SBATCH --partition=all
#SBATCH --mem-per-cpu=6000
#SBATCH --ntasks=1
#SBATCH --time=03:30:00
#SBATCH --output=/home/bandersen/iris-home/%j.out
#
#DONT FORGET
#the directory for the batch output (above) has to exist or the job won't run
cd /home/bandersen/iris-home
module load license_intel intel hdf hdf5 netcdf4 miniconda
#sh codes/dr_1_copy_class_data.sh

export HS_RET_DIR=/apollo/jung/bandersen/CSPP_UW_HSRTV_2_0
source $HS_RET_DIR/env/uw_hs_l2.bash_env
run_HSRTV.scr 3 /home/bandersen/iris-home/temp_dr_directory

python codes/dr_3_output_to_txt.py

#srun executable if it was set up for parallel processing
#otherwise just running the python like at the command line
#will have the resources of the batch job you requested
#at the top of the script




