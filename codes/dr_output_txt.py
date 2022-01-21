#Run before in terminal: module load license_intel intel hdf hdf5 netcdf4 miniconda

###	This code is to make a single txt file with all the pertinate dual regression information	###
###	for use in comparing the different methods.							###

import matplotlib.pyplot as plt
import numpy as np
import h5py
import sys
import cartopy
import cartopy.crs as ccrs
from numpy import loadtxt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from datetime import date
from shutil import copyfile
from glob import glob

def dr_txt(st1,st2,sp1,dt1,dp1,dd1,mon1,yer1):
 b1=int(st1)
 b2=int(st2)
 b3=int(sp1)
 b4=int(dt1)
 b5=int(dp1)
 b6=int(dd1)
 b7=int(mon1)
 b8=int(yer1)
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
 if b6>=0 and b6<=9:
  c6="0"+str(b6)
 else:
  c6=str(b6)
 if b7>=0 and b7<=9:
  c7="0"+str(b7)
 else:
  c7=str(b7)
 c8=str(b8)
 st=c1
 fst=c2
 sp=c3
 dt=c4
 dp=c5
 dd=c6
 mon=c7
 yer=c8
 output=('/apollo/jung/bandersen/txt_dr_output_'+yer+'_'+mon+'_'+dd+'/dr_'+yer+'_'+mon+'_'+dd+'_'+st+'_'+sp+'.txt')
 files1=glob('/apollo/jung/bandersen/raw_dr_output_'+yer+'_'+mon+'_'+dd+'/CrIS_FSR_j01_d'+yer+mon+dt+'_t'+fst+'*.h5')[:]
 files2=glob('/apollo/jung/bandersen/raw_dr_output_'+yer+'_'+mon+'_'+dd+'/CrIS_FSR_j01_d'+yer+mon+dd+'_t'+st+'*.h5')[:]
 files3=glob('/apollo/jung/bandersen/raw_dr_output_'+yer+'_'+mon+'_'+dd+'/CrIS_FSR_j01_d'+yer+mon+dp+'_t'+sp+'*.h5')[:]
 if np.size(files1)!=0 and np.size(files2)!=0 and np.size(files3)!=0:
  files=np.concatenate((files1[-1],files2[:],files3[0]))[:]
 elif np.size(files1)!=0 and np.size(files2)!=0 and np.size(files3)==0:
  files=np.concatenate((files1[-1],files2[:]))[:]
 elif np.size(files1)==0 and np.size(files2)!=0 and np.size(files3)!=0:
  files=np.concatenate((files2[:],files3[0]))[:]
 else:
  files=files2[:]
 sizer=np.size(files)
 latt=[]
 lonn=[]
 cott=[]
 ctpp=[]
 cmss=[]
 for i in range(0,sizer):
  ff=files[i]
  f1=h5py.File(ff,'r')
  shape=np.size(f1['Cmask'][:])
  cmskf1=np.reshape(f1['Cmask'],(shape))
  cmsf1=cmskf1[np.where(cmskf1!=-9999)]
  ctppf1=np.reshape(f1['CTP'][:],(shape))
  ctpf1=ctppf1[np.where(cmskf1!=-9999)]
  cottf1=np.reshape(f1['COT'][:],(shape))
  cotf1=cottf1[np.where(cmskf1!=-9999)]
  lattf1=np.reshape(f1['Latitude'][:],(shape))
  latf1=lattf1[np.where(cmskf1!=-9999)]
  lonnf1=np.reshape(f1['Longitude'][:],(shape))
  lonf1=lonnf1[np.where(cmskf1!=-9999)]
  print(np.shape(lonf1))
  latt=np.concatenate((latt,latf1))
  lonn=np.concatenate((lonn,lonf1))
  cott=np.concatenate((cott,cotf1))
  ctpp=np.concatenate((ctpp,ctpf1))
  cmss=np.concatenate((cmss,cmsf1))
  print(np.shape(lonn))
 print(ctpp)
 print(np.shape(ctpp))
 size=np.size(latt)
 lat=np.reshape(latt,(size))
 lon=np.reshape(lonn,(size))
 cot=np.reshape(cott,(size))
 ctp=np.reshape(ctpp,(size))
 cms=np.reshape(cmss,(size))
 print(ctp)
 print(np.shape(ctp))
 file=open(output,"w+")
 with open(output,'wb') as f:
  for a in range(0,size):
   z=np.array([lat[a],lon[a],cot[a],ctp[a],cms[a]])
   np.savetxt(f,z,fmt='%1.7f',newline=" ")
   np.savetxt(f,[-8888],fmt='%d',newline=" ")
 f.close()
 f=open(output,"rt")
 r=f.read()
 r=r.replace("-8888","\n")
 r=r.replace(" -","-")
 r=r.replace("ift","888888")
 r=r.replace("nan","-9999")
 r=r.replace(" ",",")
 r=r.replace(","+"\n","\n")
 r=r.replace("\n"+",","\n")
 r=r.replace("0-","0,-")
 r=r.replace("1-","1,-")
 r=r.replace("2-","2,-")
 r=r.replace("3-","3,-")
 r=r.replace("4-","4,-")
 r=r.replace("5-","5,-")
 r=r.replace("6-","6,-")
 r=r.replace("7-","7,-")
 r=r.replace("8-","8,-")
 r=r.replace("9-","9,-")
 f.close()
 f=open(output,"wt")
 f.write(r)
 f.close()
