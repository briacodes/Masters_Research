#Run before in terminal: module load license_intel intel hdf hdf5 netcdf4 miniconda

###     This code runs the full CO2 slicing technique then outputs a txt file with the information.     ###
###     It also provides the colocated VIIRS data in the output file.                                   ###

from datetime import datetime
start_time=datetime.now()
import matplotlib
matplotlib.use("AGG")
import sys
sys.path.insert(1,'/home/bandersen/iris-home/all_algorithms/functions')
import ecmwf_cd_funcs as func
import subprocess
#so for 2021-04-11 0000-0030 you would put ecmwf_full(00,00,00,30,11,11)
a=0
b=1
c=30
d=30
e=7
f=2021
func.ecmwf_full(a,b,c,d,e,f)
end_time=datetime.now()

print("done")
end_time=datetime.now()
print('Duration: {}'.format(end_time-start_time))
