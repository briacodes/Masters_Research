#Run before in terminal: module load license_intel intel hdf hdf5 netcdf4 miniconda

###     There are two forms of this code. One is to give both the categorical FOV maps of co-located    ###
###     VIIRS, CO2 Slicing, and DR as well as the straight up COT and CTP maps. This code also has      ###
###     the framework for historgram comparative imagaes.                                               ###


from datetime import datetime
start_time=datetime.now()

import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt
import numpy as np
#import h5py as h
import os
import sys
import cartopy
import cartopy.crs as ccrs
import matplotlib.colors as colors
import matplotlib.patches as patches
import matplotlib.ticker as mticker
from numpy import loadtxt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from datetime import date
from shutil import copyfile
import shutil
import netCDF4 as cdf
import sys
sys.path.insert(1,'/home/bandersen/iris-home/all_algorithms/funcs')
import functions_co2_slicing as func

f1=cdf.Dataset("/apollo/jung/snebuda/brianne/2021/cris.431.n20.20210411_173000.20210411_180000.bt.nc")
f2=cdf.Dataset("/apollo/jung/snebuda/brianne/2021/cris.431.n20.20210411_173000.20210411_180000.wf.nc")
f3=cdf.Dataset("/apollo/jung/snebuda/brianne/2021/cris.431.n20.20210411_173000.20210411_180000.state.nc")
f4=cdf.Dataset("/apollo/jung/snebuda/brianne/2021/cris.431.n20.20210411_173000.20210411_180000.te.nc")
da=f1['date'][:]
la=f1['lat'][:]
lo=f1['lon'][:]
fo=f1['solzen'][:]
fv=f1['fovn'][:]
tc=f1['tocc'][:]
to=f1['lsql'][:]
lf=f1['alfr'][:]
pv=f1['hoct'][:]
bi=f1['bias_gfs'][:]
us=f1['useflag'][:]
ob=f1['tb_obs'][:]
cr=f1['tb_clear'][:]
wf=f2['wf'][:]
pr=f3['pmid'][:]
te=f3['temp'][:]
ta=f4['ptau5'][:]


wn1=32
wn2=53
wn3=69
wn4=98
wn5=122
wn6=216
num=1585

pbl_top=pr[num,func.pbl(num,pr,te)]
tropp=pr[num,func.trop(num,pr,te,da,la)]
psp=np.size(func.integrate_six(num,pr,te,da,la)[:])-2
presures=func.integrate_six(num,pr,te,da,la)[2:psp]
gr1=func.integrate_five(num,wn1,wn2,ob,cr,ta,te,pr,da,la,bi)
gr2=func.integrate_five(num,wn2,wn3,ob,cr,ta,te,pr,da,la,bi)
gr3=func.integrate_five(num,wn3,wn4,ob,cr,ta,te,pr,da,la,bi)
gr4=func.integrate_five(num,wn4,wn5,ob,cr,ta,te,pr,da,la,bi)

a1_01=func.co2(num,28,cr,ob,bi)
a1_02=func.co2(num,30,cr,ob,bi)
a1_03=func.co2(num,32,cr,ob,bi)
a1_04=func.co2(num,34,cr,ob,bi)
a1_05=func.co2(num,36,cr,ob,bi)
a1_06=func.co2(num,49,cr,ob,bi)
a1_07=func.co2(num,51,cr,ob,bi)
a1_08=func.co2(num,53,cr,ob,bi)
a1_09=func.co2(num,55,cr,ob,bi)
a1_10=func.co2(num,57,cr,ob,bi)
a1_11=func.co2(num,65,cr,ob,bi)
a1_12=func.co2(num,67,cr,ob,bi)
a1_13=func.co2(num,69,cr,ob,bi)
a1_14=func.co2(num,71,cr,ob,bi)
a1_15=func.co2(num,73,cr,ob,bi)
a1_16=func.co2(num,94,cr,ob,bi)
a1_17=func.co2(num,96,cr,ob,bi)
a1_18=func.co2(num,98,cr,ob,bi)
a1_19=func.co2(num,100,cr,ob,bi)
a1_20=func.co2(num,102,cr,ob,bi)
a1_21=func.co2(num,122,cr,ob,bi)

print(a1_03)
print(a1_08)
print(a1_13)
print(a1_18)

fig=plt.figure()
plt.plot(gr4,presures,color='purple',label=r"$A_3(733.125,748.125)$")
plt.axhline(y=tropp,color='black',linestyle=":",label="Tropopause Height")
plt.axvline(x=0,color='black')
plt.axhline(y=pbl_top,color='black',linestyle="-.",label="PBL Top Height")
plt.ylabel("Pressure hPa")
plt.xlim(-0.5,0.5)
plt.yscale('log',basey=2)
plt.yticks([150,200,300,400,500,600,700,800,900,1000],['150','200','300','400','500','600','700','800','900','1000'])
plt.gca().invert_yaxis()
plt.xlabel(r"$A_{3}$ Value at Pressure Level")
plt.ylim(1000,150)
filepth=("/home/bandersen/iris-home/co2_exampl_plot_2.png")
fig.set_size_inches(10,10)
plt.savefig(filepth)

xs1=[a1_03]
xs2=[a1_08]
xs3=[a1_13]
xs4=[a1_18]
xs5=[a1_21]
ys=[0.5]

fig=plt.figure()
plt.scatter(xs1,ys,color='red',label=r"$A_1(691.875)$")
plt.scatter(xs2,ys,color='blue',label=r"$A_1(705)$")
plt.scatter(xs3,ys,color='green',label=r"$A_1(715)$")
plt.scatter(xs4,ys,color='purple',label=r"$A_1(733.125)$")
plt.scatter(xs5,ys,color='grey',label=r"$A_1(748.125)$")
plt.axvline(x=0.5,color='black',label="Cutoff")
plt.ylim(0,1)
plt.xlim(-0.4,1.4)
filepth=("/home/bandersen/iris-home/co2_exampl_plot_1.png")
fig.set_size_inches(10,2)
plt.savefig(filepth)


