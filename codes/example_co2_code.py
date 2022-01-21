#Run before in terminal: module load license_intel intel hdf hdf5 netcdf4 miniconda

###     This code runs the full CO2 slicing technique then outputs a txt file with the information.     ###
###     It also provides the colocated VIIRS data in the output file.                                   ###

from datetime import datetime
start_time=datetime.now()
import matplotlib
matplotlib.use("AGG")
import sys
sys.path.insert(1,'/home/bandersen/iris-home/all_algorithms/functions')
import comp_func as comp_func
import functions_co2_slicing_h_diff as cs_func
#so for 2021-04-11 0000-0030 you would put ecmwf_full(00,00,00,30,11,11)
a=0
b=1
c=10
d=10
e=10
f=2020

cs_func.full_thing(a,0,a,30,c,c,e,f,0.5)
cs_func.full_thing(a,30,b,0,c,d,e,f,0.5)
comp_func.full_comps_parts(a,b,c,d,e,f)

print("done")
end_time=datetime.now()
print('Duration: {}'.format(end_time-start_time))


