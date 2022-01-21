#Run before in terminal: module load license_intel intel hdf hdf5 netcdf4 miniconda

###	This code it is to take netcdf files with given cloud info and transform it into a formate which	###
###	can be used in the ECMWF Cloud Detect Software using McNally & Watts work.				###

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
sys.path.insert(1,'/home/bandersen/iris-home/full_co2_slicing_2021_04/cs')
import functions_co2_slicing_h_diff as cs_func
today=date.today()

#so for 2021-04-10 0000-0030 you would put ecmwf_full(00,00,00,30,10,10)
def ecmwf_full(st1,sp1,dt1,dp1,mon1,yer1):
 f1=cdf.Dataset("/apollo/jung/snebuda/brianne/2021/cris.431.n20.20210410_000000.20210410_003000.bt.nc")
 chan=f1['channum'][:]
 wave=[]
 for i in range(0,431):
  z=chan[i]
  if z<=713:
   a=650+(.625*(z-1))
   wave.append(a)
  elif z>=714 and z<=1578:
   a=1210+(.625*(z-714))
   wave.append(a)
  elif z>=1579:
   a=2155+(.625*(z-1579))
   wave.append(a)
  else:
   pass
 chan_size=np.size(chan)
 b1=int(st1)
 b2=int(sp1)
 b3=int(dt1)
 b4=int(dp1)
 b5=int(mon1)
 b6=int(yer1)
 if b1>=0 and b1<=9:
  c1="0"+str(b1)
 else:
  c1=str(b1)
 if b2>=0 and b2<=9:
  c2="0"+str(b2)
 else:
  c2=str(b2)
 if b3>=0 and b3<=9:
  c3="0"+str(b3)
 else:
  c3=str(b3)
 if b4>=0 and b4<=9:
  c4="0"+str(b4)
 else:
  c4=str(b4)
 if b5>=0 and b5<=9:
  c5="0"+str(b5)
 else:
  c5=str(b5)
 c6=str(b6)
 st=c1
 sp=c2
 dt=c3
 dp=c4
 mon=c5
 yer=c6
 fa1=cdf.Dataset("/apollo/jung/snebuda/brianne/2021/cris.431.n20."+yer+mon+dt+"_"+st+"0000."+yer+mon+dt+"_"+st+"3000.bt.nc")
 fa2=cdf.Dataset("/apollo/jung/snebuda/brianne/2021/cris.431.n20."+yer+mon+dt+"_"+st+"0000."+yer+mon+dt+"_"+st+"3000.wf.nc")
 fa3=cdf.Dataset("/apollo/jung/snebuda/brianne/2021/cris.431.n20."+yer+mon+dt+"_"+st+"0000."+yer+mon+dt+"_"+st+"3000.state.nc")
 fb1=cdf.Dataset("/apollo/jung/snebuda/brianne/2021/cris.431.n20."+yer+mon+dt+"_"+st+"3000."+yer+mon+dp+"_"+sp+"0000.bt.nc")
 fb2=cdf.Dataset("/apollo/jung/snebuda/brianne/2021/cris.431.n20."+yer+mon+dt+"_"+st+"3000."+yer+mon+dp+"_"+sp+"0000.wf.nc")
 fb3=cdf.Dataset("/apollo/jung/snebuda/brianne/2021/cris.431.n20."+yer+mon+dt+"_"+st+"3000."+yer+mon+dp+"_"+sp+"0000.state.nc")
 day=np.concatenate((fa1['date'][1:],fb1['date'][1:]),axis=0)
 lat=np.concatenate((fa1['lat'][1:],fb1['lat'][1:]),axis=0)
 lon=np.concatenate((fa1['lon'][1:],fb1['lon'][1:]),axis=0)
 pres=np.concatenate((fa3['pmid'][1:],fb3['pmid'][1:]),axis=0)
 tb_clear_a=np.add(fa1['tb_clear'][1:],fa1['bias_gfs'][:])
 tb_clear_b=np.add(fb1['tb_clear'][1:],fb1['bias_gfs'][:])
 tb_clear=np.concatenate((tb_clear_a,tb_clear_b),axis=0)
 tb_obs=np.concatenate((fa1['tb_obs'][1:],fb1['tb_obs'][1:]),axis=0)
 temp=np.concatenate((fa3['temp'][1:],fb3['temp'][1:]),axis=0)
 wf=np.concatenate((fa2['wf'][1:],fb2['wf'][1:]),axis=0)
 fov_size=np.size(lat)
 chan_size=np.size(chan)
 outfile=("/apollo/jung/bandersen/ecmwf_data_input_"+yer+"_"+mon+"_"+dt+"/ecmwf_cd_input_"+yer+mon+dt+"_n20_"+st+"00-"+sp+"00.txt")
 heights=[]
 for ii in range(0,fov_size):
  for i in range(0,chan_size):
   heights.append(pres[ii][np.where(wf[ii][i]==np.max(wf[ii][i]))[0]][0])
 heights=np.reshape(heights,(fov_size,chan_size))
 with open(outfile,'wb') as f:
  np.savetxt(f,[27],fmt='%d',newline=" ")
  np.savetxt(f,[-8888],fmt='%d',newline=" ")
  np.savetxt(f,[int(chan_size)],fmt='%d',newline=" ")
  np.savetxt(f,[-8888],fmt='%d',newline=" ")
  np.savetxt(f,chan,fmt='%d',newline=" ")
  np.savetxt(f,[-8888],fmt='%d',newline=" ")
  np.savetxt(f,[int(fov_size)],fmt='%d',newline=" ")
  for i in range(0,fov_size):
   a=np.array([lon[i],lat[i]])
   b=np.array([pres[cs_func.trop(i,pres,temp,day,lat)],pres[cs_func.pbl(i,pres,temp)]])
   c=tb_obs[i][:]
   d=tb_clear[i][:]
   e=heights[i][:]
   np.savetxt(f,[-8888],fmt='%d',newline=" ")
   np.savetxt(f,a,fmt='%1.5f',newline=" ")
   np.savetxt(f,[-8888],fmt='%d',newline=" ")
   np.savetxt(f,b,fmt='%d',newline=" ")
   np.savetxt(f,[-8888],fmt='%d',newline=" ")
   np.savetxt(f,c,fmt='%1.5f',newline=" ")
   np.savetxt(f,[-8888],fmt='%d',newline=" ")
   np.savetxt(f,d,fmt='%1.5f',newline=" ")
   np.savetxt(f,[-8888],fmt='%d',newline=" ")
   np.savetxt(f,e,fmt='%1.5f',newline=" ")
 f.close()
 f=open(outfile,"rt")
 r=f.read()
 r=r.replace("-8888","\n")
 f.close()
 f=open(outfile,"wt")
 f.write(r)
 f.close()
 return("done")

def file_name(st1,sp1,dt1,mon1,yer1):
 b1=int(st1)
 b2=int(sp1)
 b3=int(dt1)
 b4=int(mon1)
 b5=int(yer1)
 if b1>=0 and b1<=9:
  c1="0"+str(b1)
 else:
  c1=str(b1)
 if b2>=0 and b2<=9:
  c2="0"+str(b2)
 else:
  c2=str(b2)
 if b3>=0 and b3<=9:
  c3="0"+str(b3)
 else:
  c3=str(b3)
 if b4>=0 and b4<=9:
  c4="0"+str(b4)
 else:
  c4=str(b4)
 c5=str(b5)
 st=c1
 sp=c2
 dt=c3
 mon=c4
 yer=c5
 filer=("/home/bandersen/iris-home/all_algorithms/"+yer+"_"+mon+"_"+dt+"_ecmwf/"+st+"_"+sp+".sh")
 fuller=(filer)
 return(fuller)

def last(st1,sp1,dt1,mon1,yer1):
 b1=int(st1)
 b2=int(sp1)
 b3=int(dt1)
 b4=int(mon1)
 b5=int(yer1)
 if b1>=0 and b1<=9:
  c1="0"+str(b1)
 else:
  c1=str(b1)
 if b2>=0 and b2<=9:
  c2="0"+str(b2)
 else:
  c2=str(b2)
 if b3>=0 and b3<=9:
  c3="0"+str(b3)
 else:
  c3=str(b3)
 if b4>=0 and b4<=9:
  c4="0"+str(b4)
 else:
  c4=str(b4)
 c5=str(b5)
 st=c1
 sp=c2
 dt=c3
 mon=c4
 yer=c5
 original_1=('/apollo/jung/bandersen/ecmwf_data_output_raw_'+yer+'_'+mon+'_'+dt+'/ecmwf_cd_output_'+yer+mon+dt+'_'+st+'00_'+sp+'00_n20.txt.newest')
 file2=('/apollo/jung/bandersen/ecmwf_data_output_txt_'+yer+'_'+mon+'_'+dt+'/ecmwf_cd_'+yer+mon+dt+'_'+st+'_'+sp+'_n20.txt')
 filenames=[original_1]
 with open(file2,'w') as outfile:
  for fname in filenames:
   with open(fname) as infile:
    for line in infile:
     outfile.write(line)
 f=open(file2,"rt")
 r=f.read()
 r=r.replace("-----","-8888")
 r=r.replace("\n"," ")
 r=r.replace("Longitude Latitude ObNumber"," ")
 r=r.replace("Cloud flags"," ")
 r=r.replace("Aerosol flags"," ")
 r=r.replace("       "," ")
 r=r.replace("      "," ")
 r=r.replace("     "," ")
 r=r.replace("    "," ")
 r=r.replace("   "," ")
 r=r.replace("  "," ")
 r=r.replace(" ",",")
 r=r.replace(",-8888,","\n")
 r=r.replace("-8888,","")
 r=r.replace(","+"\n","\n")
 f.close()
 f=open(file2,"wt")
 f.write(r)
 f.close()
 f=open(file2,"a")
 f.write("-7777")
 f.close()
 f=open(file2,"rt")
 r=f.read()
 r=r.replace(",-7777"," ")
 f.close()
 f=open(file2,"wt")
 f.write(r)
 f.close()
 return("thi")
