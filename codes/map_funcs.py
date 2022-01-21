#Run before in terminal: module load license_intel intel hdf hdf5 netcdf4 miniconda

###     This code provided the functions used to make the map type plots shown in the thesis and 	###
###     analysis for Brianne Andersen's Master's work. 							###

from datetime import datetime
start_time=datetime.now()

import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt
import numpy as np
import h5py as h
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

extent=[-180,180,-90,90]
x_ticks1=[-180,-165,-150,-135,-120,-105,-90,-75,-60,-45,-30,-15,0,15,30,45,60,75,90,105,120,135,150,165,180]
y_ticks1=[-90,-75,-60,-45,-30,-15,0,15,30,45,60,75,90]
x_ticks2=[-180,-135,-90,-45,0,45,90,135,180]
y_ticks2=[-90,-60,-30,0,30,60,90]

#Type codes are as follows:
#1=CO2 Slicing maps
#2=NCEP Algorithm
#3=DR
#this code can be modified to exclude or include certain boundaries (such as CTP above 300 hPa) by adding addtional
#np.where(ctpc<=50,2500,ctpc) phrases
def ctp_maps(yr,mn,dy,type):
 b1=int(dy)
 b2=int(mn)
 b3=int(yr)
 if b1>=0 and b1<=9:
  c1="0"+str(b1)
 else:
  c1=str(b1)
 if b2>=0 and b2<=9:
  c2="0"+str(b2)
 else:
  c2=str(b2)
 c3=str(b3)
 day=c1
 mon=c2
 yer=c3
 day_time=(yer+'_'+mon+'_'+day)
 if type==1:
  file=("/apollo/jung/bandersen/compare_full_no_parts/co2_ncep_output_"+day_time+"_24hr.txt")
  array=loadtxt(file,comments='#',delimiter=',',unpack=True)
  filepath=("/home/bandersen/iris-home/all_algorithms/image_output")
  day_label=(r+day_time+'\n')
  lat=array[0][:]
  lon=array[1][:]
  ctpc=array[9][:]
  ctpc=np.where(ctpc<=50,2500,ctpc)
  ctpc[ctpc>=850]=np.nan
  ctpc_cld=[]
  for i in range(0,size):
   a=ctpc[i]
   ctpc_cld.append(a)
  ctpc=array[9][:]
  ctpc=np.where(ctpc<=1900,2200,ctpn)
  ctpc[ctpc>=2050]=np.nan
  ctpc_no_cld_co2=[]
  for i in range(0,size):
   a=ctpc[i]
   ctpc_no_cld_co2.append(a)
  ctpc=array[9][:]
  ctpc=np.where(ctpc<=1999,0,ctpc)
  ctpc[ctpc<=0]=np.nan
  ctpc_error=[]
  for i in range(0,size):
   a=ctpc[i]
   ctpc_error.append(a)
  fig=plt.figure()
  proj=ccrs.PlateCarree()
  ax=plt.axes(projection=proj)
  ax.add_feature(cartopy.feature.OCEAN,zorder=0,alpha=0.5,facecolor='lightskyblue')
  ax.add_feature(cartopy.feature.BORDERS,edgecolor='slategray')
  ax.add_feature(cartopy.feature.COASTLINE,edgecolor='black')
  ax.add_feature(cartopy.feature.LAND,zorder=0,edgecolor='black',alpha=0.5,facecolor='beige')
  ax.set_extent(extent,crs=proj)
  gl=ax.gridlines(crs=proj,linewidth=2,color='grey',alpha=0.5,linestyle='--',draw_labels=True)
  gl.xlocator=mticker.FixedLocator(x_ticks1)
  gl.ylocator=mticker.FixedLocator(y_ticks1)
  gl.xlabels_top=False
  gl.xlabels_bottom=False
  gl.ylabels_right=False
  gl.ylabels_left=False
  gl.xlines=True
  gl.xformatter=LONGITUDE_FORMATTER
  gl.yformatter=LATITUDE_FORMATTER
  gl.xlabel_style={'color':'red','weight':'bold','size':'40'}
  gl.ylabel_style={'color':'red','weight':'bold','size':'40'}
  g2=ax.gridlines(crs=proj,linewidth=2,color='grey',alpha=0.5,linestyle='--',draw_labels=True)
  g2.xlocator=mticker.FixedLocator(x_ticks2)
  g2.ylocator=mticker.FixedLocator(y_ticks2)
  g2.xlabels_bottom=False
  g2.ylabels_right=False
  g2.ylabels_left=True
  g2.xlines=True
  g2.xformatter=LONGITUDE_FORMATTER
  g2.yformatter=LATITUDE_FORMATTER
  g2.xlabel_style={'color':'black','weight':'bold','size':'14'}
  g2.ylabel_style={'color':'black','weight':'bold','size':'14'}
  plt.scatter(lon,lat,c=ctpc_no_cld_co2,cmap='gist_yarg',s=0.05,vmin=1999,vmax=5000)
  data=plt.scatter(lon,lat,c=ctpc_cld,cmap='gist_rainbow',s=0.05,vmin=50,vmax=950)
  plt.title('24-hour CO2 Slicing Cloud Top Pressure Calculation'+'\n'+day_label+'\n',fontsize=22,weight='bold')
  cb=plt.colorbar(data,aspect=25,fraction=.025,pad=.04,orientation='horizontal')
  cb.set_label(label='\nPressure [hPa]',fontsize=18)
  cb.ax.tick_params(labelsize=16)
  filepth=(filepath+"/ctp_co2_"+day_time+".png")
  fig.set_size_inches(12,8.75)
  plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight')
 elif type==2:
  file=("/apollo/jung/bandersen/compare_full_no_parts/co2_ncep_output_"+day_time+"_24hr.txt")
  array=loadtxt(file,comments='#',delimiter=',',unpack=True)
  filepath=("/home/bandersen/iris-home/all_algorithms/image_output")
  day_label=(r+day_time+'\n')
  lat=array[0][:]
  lon=array[1][:]
  ctpn=array[11][:]
  ctpn=np.where(ctpn<=50,100,ctpn)
  ctpn[ctpn>=850]=np.nan
  ctpn_cld=[]
  for i in range(0,size):
   a=ctpn[i]
   ctpn_cld.append(a)
  ctpn=array[11][:]
  ctpn=np.where(ctpn<=1900,2200,ctpn)
  ctpn[ctpn>=2050]=np.nan
  ctpn_no_cld_ncep=[]
  for i in range(0,size):
   a=ctpn[i]
   ctpn_no_cld_ncep.append(a)
  ctpn=array[11][:]
  ctpn=np.where(ctpn<=1999,0,ctpn)
  ctpn[ctpn<=0]=np.nan
  ctpn_error=[]
  for i in range(0,size):
   a=ctpn[i]
   ctpn_error.append(a)
  fig=plt.figure()
  proj=ccrs.PlateCarree()
  ax=plt.axes(projection=proj)
  ax.add_feature(cartopy.feature.OCEAN,zorder=0,alpha=0.5,facecolor='lightskyblue')
  ax.add_feature(cartopy.feature.BORDERS,edgecolor='slategray')
  ax.add_feature(cartopy.feature.COASTLINE,edgecolor='black')
  ax.add_feature(cartopy.feature.LAND,zorder=0,edgecolor='black',alpha=0.5,facecolor='beige')
  ax.set_extent(extent,crs=proj)
  gl=ax.gridlines(crs=proj,linewidth=2,color='grey',alpha=0.5,linestyle='--',draw_labels=True)
  gl.xlocator=mticker.FixedLocator(x_ticks1)
  gl.ylocator=mticker.FixedLocator(y_ticks1)
  gl.xlabels_top=False
  gl.xlabels_bottom=False
  gl.ylabels_right=False
  gl.ylabels_left=False
  gl.xlines=True
  gl.xformatter=LONGITUDE_FORMATTER
  gl.yformatter=LATITUDE_FORMATTER
  gl.xlabel_style={'color':'red','weight':'bold','size':'40'}
  gl.ylabel_style={'color':'red','weight':'bold','size':'40'}
  g2=ax.gridlines(crs=proj,linewidth=2,color='grey',alpha=0.5,linestyle='--',draw_labels=True)
  g2.xlocator=mticker.FixedLocator(x_ticks2)
  g2.ylocator=mticker.FixedLocator(y_ticks2)
  g2.xlabels_bottom=False
  g2.ylabels_right=False
  g2.ylabels_left=True
  g2.xlines=True
  g2.xformatter=LONGITUDE_FORMATTER
  g2.yformatter=LATITUDE_FORMATTER
  g2.xlabel_style={'color':'black','weight':'bold','size':'14'}
  g2.ylabel_style={'color':'black','weight':'bold','size':'14'}
  plt.scatter(lon,lat,c=ctpn_no_cld_ncep,cmap='gist_yarg',s=0.05,vmin=1999,vmax=5000)
  data=plt.scatter(lon,lat,c=ctpn_cld,cmap='gist_rainbow',s=0.05,vmin=50,vmax=950)
  plt.title('24-hour NCEP Algorithm Cloud Top Pressure Calculation'+'\n'+day_label+'\n',fontsize=22,weight='bold')
  cb=plt.colorbar(data,aspect=25,fraction=.025,pad=.04,orientation='horizontal')
  cb.set_label(label='\nPressure [hPa]',fontsize=18)
  cb.ax.tick_params(labelsize=16)
  filepth=(filepath+"/ctp_ncep_"+day_time+".png")
  fig.set_size_inches(12,8.75)
  plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight')
 elif type==3:
  file=("/apollo/jung/bandersen/compare_full_no_parts/co2_dr_ncep_ecmwf_output_"+day_time+"_24hr.txt")
  array=loadtxt(file,comments='#',delimiter=',',unpack=True)
  filepath=("/home/bandersen/iris-home/all_algorithms/image_output")
  day_label=(r+day_time+'\n')
  lat=array[0][:]
  lon=array[1][:]
  ctpd=array[10][:]
  ctpd=np.where(ctpd<=50,100,ctpd)
  ctpd[ctpd>=850]=np.nan
  ctpd_cld=[]
  for i in range(0,size):
   a=ctpd[i]
   ctpd_cld.append(a)
  ctpd=array[10][:]
  ctpd=np.where(ctpd<=1900,2200,ctpd)
  ctpd[ctpd>=2050]=np.nan
  ctpd_no_cld_dr=[]
  for i in range(0,size):
   a=ctpd[i]
   ctpd_no_cld_co2.append(a)
  ctpd=array[10][:]
  ctpd=np.where(ctpd<=1999,0,ctpd)
  ctpd[ctpd<=0]=np.nan
  ctpd_error=[]
  for i in range(0,size):
   a=ctpd[i]
   ctpd_error.append(a)
  fig=plt.figure()
  proj=ccrs.PlateCarree()
  ax=plt.axes(projection=proj)
  ax.add_feature(cartopy.feature.OCEAN,zorder=0,alpha=0.5,facecolor='lightskyblue')
  ax.add_feature(cartopy.feature.BORDERS,edgecolor='slategray')
  ax.add_feature(cartopy.feature.COASTLINE,edgecolor='black')
  ax.add_feature(cartopy.feature.LAND,zorder=0,edgecolor='black',alpha=0.5,facecolor='beige')
  ax.set_extent(extent,crs=proj)
  gl=ax.gridlines(crs=proj,linewidth=2,color='grey',alpha=0.5,linestyle='--',draw_labels=True)
  gl.xlocator=mticker.FixedLocator(x_ticks1)
  gl.ylocator=mticker.FixedLocator(y_ticks1)
  gl.xlabels_top=False
  gl.xlabels_bottom=False
  gl.ylabels_right=False
  gl.ylabels_left=False
  gl.xlines=True
  gl.xformatter=LONGITUDE_FORMATTER
  gl.yformatter=LATITUDE_FORMATTER
  gl.xlabel_style={'color':'red','weight':'bold','size':'40'}
  gl.ylabel_style={'color':'red','weight':'bold','size':'40'}
  g2=ax.gridlines(crs=proj,linewidth=2,color='grey',alpha=0.5,linestyle='--',draw_labels=True)
  g2.xlocator=mticker.FixedLocator(x_ticks2)
  g2.ylocator=mticker.FixedLocator(y_ticks2)
  g2.xlabels_bottom=False
  g2.ylabels_right=False
  g2.ylabels_left=True
  g2.xlines=True
  g2.xformatter=LONGITUDE_FORMATTER
  g2.yformatter=LATITUDE_FORMATTER
  g2.xlabel_style={'color':'black','weight':'bold','size':'14'}
  g2.ylabel_style={'color':'black','weight':'bold','size':'14'}
  plt.scatter(lon,lat,c=ctpd_no_cld_dr,cmap='gist_yarg',s=0.05,vmin=1999,vmax=5000)
  data=plt.scatter(lon,lat,c=ctpd_cld,cmap='gist_rainbow',s=0.05,vmin=50,vmax=950)
  plt.title('24-hour DR Cloud Top Pressure Calculation'+'\n'+day_label+'\n',fontsize=22,weight='bold')
  cb=plt.colorbar(data,aspect=25,fraction=.025,pad=.04,orientation='horizontal')
  cb.set_label(label='\nPressure [hPa]',fontsize=18)
  cb.ax.tick_params(labelsize=16)
  filepth=(filepath+"/ctp_dr_"+day_time+".png")
  fig.set_size_inches(12,8.75)
  plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight')
 else:
  p=2

#type 1 is for including dr and ecmwf files, type 2 is only ncep and co2 slicing
def ctp_ncep_co2_slicing_diff_plots(yr,mn,dy,type):
 b1=int(dy)
 b2=int(mn)
 b3=int(yr)
 if b1>=0 and b1<=9:
  c1="0"+str(b1)
 else:
  c1=str(b1)
 if b2>=0 and b2<=9:
  c2="0"+str(b2)
 else:
  c2=str(b2)
 c3=str(b3)
 day=c1
 mon=c2
 yer=c3
 day_time=(yer+'_'+mon+'_'+day)
 if type==1:
  file=("/apollo/jung/bandersen/compare_full_no_parts/co2_dr_ncep_ecmwf_output_"+day_time+"_24hr.txt")
  array=loadtxt(file,comments='#',delimiter=',',unpack=True)
  filepath=("/home/bandersen/iris-home/all_algorithms/image_output")
  day_label=(r+day_time+'\n')
  lat=array[0][:]
  lon=array[1][:]
  ctpc=array[9][:]
  ctpn=array[12][:]
  ctp=ctpc-ctpn
  ctp=np.where(ctp>=1000,-1000,ctp)
  ctp[ctp<=-1000]=np.nan
  fig=plt.figure()
  plt.hist(ctp[~np.isnan(ctp)],bins=201,range=(-1000,1000))
  plt.title("CO2 Slicing and NCEP CTP Difference [CO2 - NCEP] Distribution"+'\n')
  plt.xlabel("CTP Difference [hPa]")
  plt.ylabel("Count")
  filepth=(filepath+"/ctp_co2_ncep_diff_bar_"+day_time+".png")
  fig.set_size_inches(12,8.75)
  plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight')
  fig=plt.figure()
  proj=ccrs.PlateCarree()
  ax=plt.axes(projection=proj)
  ax.add_feature(cartopy.feature.OCEAN,zorder=0,alpha=0.5,facecolor='lightskyblue')
  ax.add_feature(cartopy.feature.BORDERS,edgecolor='slategray')
  ax.add_feature(cartopy.feature.COASTLINE,edgecolor='black')
  ax.add_feature(cartopy.feature.LAND,zorder=0,edgecolor='black',alpha=0.5,facecolor='beige')
  ax.set_extent(extent,crs=proj)
  gl=ax.gridlines(crs=proj,linewidth=2,color='grey',alpha=0.5,linestyle='--',draw_labels=True)
  gl.xlocator=mticker.FixedLocator(x_ticks1)
  gl.ylocator=mticker.FixedLocator(y_ticks1)
  gl.xlabels_top=False
  gl.xlabels_bottom=False
  gl.ylabels_right=False
  gl.ylabels_left=False
  gl.xlines=True
  gl.xformatter=LONGITUDE_FORMATTER
  gl.yformatter=LATITUDE_FORMATTER
  gl.xlabel_style={'color':'red','weight':'bold','size':'40'}
  gl.ylabel_style={'color':'red','weight':'bold','size':'40'}
  g2=ax.gridlines(crs=proj,linewidth=2,color='grey',alpha=0.5,linestyle='--',draw_labels=True)
  g2.xlocator=mticker.FixedLocator(x_ticks2)
  g2.ylocator=mticker.FixedLocator(y_ticks2)
  g2.xlabels_bottom=False
  g2.ylabels_right=False
  g2.ylabels_left=True
  g2.xlines=True
  g2.xformatter=LONGITUDE_FORMATTER
  g2.yformatter=LATITUDE_FORMATTER
  g2.xlabel_style={'color':'black','weight':'bold','size':'14'}
  g2.ylabel_style={'color':'black','weight':'bold','size':'14'}
  data=plt.scatter(lon[~np.isnan(ctp)],lat[~np.isnan(ctp)],c=ctp[~np.isnan(ctp)],cmap='seismic',s=0.05,vmin=-1000,vmax=1000)
  plt.title('24-hour CO2 Slicing Minus NCEP Algorithm'+'\n'+'Cloud Top Pressure Calculation '+day_label+'\n',fontsize=22,weight='bold')
  cb=plt.colorbar(data,aspect=25,fraction=.025,pad=.04,orientation='horizontal')
  cb.set_label(label='\nPressure [hPa]',fontsize=18)
  cb.ax.tick_params(labelsize=16)
  filepth=(filepath+"/ctp_co2_ncep_diff_"+day_time+".png")
  fig.set_size_inches(12,8.75)
  plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight')
 elif type==2:
  file=("/apollo/jung/bandersen/compare_full_no_parts/co2_ncep_output_"+day_time+"_24hr.txt")
  array=loadtxt(file,comments='#',delimiter=',',unpack=True)
  filepath=("/home/bandersen/iris-home/all_algorithms/image_output")
  day_label=(r+day_time+'\n')
  lat=array[0][:]
  lon=array[1][:]
  ctpc=array[9][:]
  ctpn=array[12][:]
  ctp=ctpc-ctpn
  ctp=np.where(ctp>=1000,-1000,ctp)
  ctp[ctp<=-1000]=np.nan
  fig=plt.figure()
  plt.hist(ctp[~np.isnan(ctp)],bins=201,range=(-1000,1000))
  plt.title("CO2 Slicing and NCEP CTP Difference [CO2 - NCEP] Distribution"+'\n')
  plt.xlabel("CTP Difference [hPa]")
  plt.ylabel("Count")
  filepth=(filepath+"/ctp_co2_ncep_diff_bar_"+day_time+".png")
  fig.set_size_inches(12,8.75)
  plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight')
  fig=plt.figure()
  proj=ccrs.PlateCarree()
  ax=plt.axes(projection=proj)
  ax.add_feature(cartopy.feature.OCEAN,zorder=0,alpha=0.5,facecolor='lightskyblue')
  ax.add_feature(cartopy.feature.BORDERS,edgecolor='slategray')
  ax.add_feature(cartopy.feature.COASTLINE,edgecolor='black')
  ax.add_feature(cartopy.feature.LAND,zorder=0,edgecolor='black',alpha=0.5,facecolor='beige')
  ax.set_extent(extent,crs=proj)
  gl=ax.gridlines(crs=proj,linewidth=2,color='grey',alpha=0.5,linestyle='--',draw_labels=True)
  gl.xlocator=mticker.FixedLocator(x_ticks1)
  gl.ylocator=mticker.FixedLocator(y_ticks1)
  gl.xlabels_top=False
  gl.xlabels_bottom=False
  gl.ylabels_right=False
  gl.ylabels_left=False
  gl.xlines=True
  gl.xformatter=LONGITUDE_FORMATTER
  gl.yformatter=LATITUDE_FORMATTER
  gl.xlabel_style={'color':'red','weight':'bold','size':'40'}
  gl.ylabel_style={'color':'red','weight':'bold','size':'40'}
  g2=ax.gridlines(crs=proj,linewidth=2,color='grey',alpha=0.5,linestyle='--',draw_labels=True)
  g2.xlocator=mticker.FixedLocator(x_ticks2)
  g2.ylocator=mticker.FixedLocator(y_ticks2)
  g2.xlabels_bottom=False
  g2.ylabels_right=False
  g2.ylabels_left=True
  g2.xlines=True
  g2.xformatter=LONGITUDE_FORMATTER
  g2.yformatter=LATITUDE_FORMATTER
  g2.xlabel_style={'color':'black','weight':'bold','size':'14'}
  g2.ylabel_style={'color':'black','weight':'bold','size':'14'}
  data=plt.scatter(lon[~np.isnan(ctp)],lat[~np.isnan(ctp)],c=ctp[~np.isnan(ctp)],cmap='seismic',s=0.05,vmin=-1000,vmax=1000)
  plt.title('24-hour CO2 Slicing Minus NCEP Algorithm'+'\n'+'Cloud Top Pressure Calculation '+day_label+'\n',fontsize=22,weight='bold')
  cb=plt.colorbar(data,aspect=25,fraction=.025,pad=.04,orientation='horizontal')
  cb.set_label(label='\nPressure [hPa]',fontsize=18)
  cb.ax.tick_params(labelsize=16)
  filepth=(filepath+"/ctp_co2_ncep_diff_"+day_time+".png")
  fig.set_size_inches(12,8.75)
  plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight')
 else:
  p=2

#for the confusion matrix plots with input of which file type you have
#type 1 is for including dr and ecmwf files, type 2 is only ncep and co2 slicing
#also with a wavenumber specific qualifier of
#wn=1=691.875
#wn=2=705
#wn=3=715
#wn=4=733.125
#wn=5=748.125
#wn=6=959.375
#green  c=1     co2     clear   ncep    clear
#blue   c=2     co2     cloud   ncep    cloud
#red    c=3     co2     clear   ncep    cloud
#purple c=4     co2     cloud   ncep    clear
def confusion_map_plots(yr,mn,dy,type,wn):
 b1=int(dy)
 b2=int(mn)
 b3=int(yr)
 if b1>=0 and b1<=9:
  c1="0"+str(b1)
 else:
  c1=str(b1)
 if b2>=0 and b2<=9:
  c2="0"+str(b2)
 else:
  c2=str(b2)
 c3=str(b3)
 day=c1
 mon=c2
 yer=c3
 day_time=(yer+'_'+mon+'_'+day)
 if type==1:
  file=("/apollo/jung/bandersen/compare_full_no_parts/co2_dr_ncep_ecmwf_output_"+day_time+"_24hr.txt")
  array=loadtxt(file,comments='#',delimiter=',',unpack=True)
  filepath=("/home/bandersen/iris-home/all_algorithms/image_output")
  day_label=(r+day_time+'\n')
  lat=array[0][:]
  lon=array[1][:]
  size=np.size(lat)
  if wn==1:
   co2=array[19][:]
   ncep=array[40][:]
   waver="691.875"
   coder="32"
  elif wn==2:
   co2=array[24][:]
   ncep=array[45][:]
   waver="705"
   coder="53"
  elif wn==3:
   co2=array[29][:]
   ncep=array[50][:]
   waver="715"
   coder="69"
  elif wn==4:
   co2=array[34][:]
   ncep=array[55][:]
   waver="733.125"
   coder="98"
  elif wn==5:
   co2=array[37][:]
   ncep=array[58][:]
   waver="748.125"
   coder="122"
  elif wn==6:
   co2=array[9][:]
   ncep=array[12][:]
   waver="Full"
   coder="full"
  else:
   print(error)
  lon1=[]
  lat1=[]
  lon2=[]
  lat2=[]
  lon3=[]
  lat3=[]
  lon4=[]
  lat4=[]
  if wn>=5:
   co2=np.where(co2==1,6,co2)
   co2=np.where(co2<=4,5,co2)
   for i in range(0,size):
    a=co2[i]
    b=ncep[i]
    if (a==6) and (b==0):
     c=lon[i]
     d=lat[i]
     lon1.append(c)
     lat1.append(d)
    elif (a==5) and (b==7):
     c=lon[i]
     d=lat[i]
     lon2.append(c)
     lat2.append(d)
    elif (a==6) and (b==7):
     c=lon[i]
     d=lat[i]
     lon3.append(c)
     lat3.append(d)
    elif (a==5) and (b==0):
     c=lon[i]
     d=lat[i]
     lon4.append(c)
     lat4.append(d)
    else:
     c=0
     d=0
  elif wn==6:
   co2=np.where(co2==2050,2000,co2)
   co2=np.where(co2==2100,300,co2)
   for i in range(0,size):
    a=co2[i]
    b=ncep[i]
    if (a>=900) and (b>=900):
     c=lon[i]
     d=lat[i]
     lon1.append(c)
     lat1.append(d)
    elif (a<900) and (b<900):
     c=lon[i]
     d=lat[i]
     lon2.append(c)
     lat2.append(d)
    elif (a>=900) and (b<900):
     c=lon[i]
     d=lat[i]
     lon3.append(c)
     lat3.append(d)
    elif (a<900) and (b>=900):
     c=lon[i]
     d=lat[i]
     lon4.append(c)
     lat4.append(d)
    else:
     c=0
     d=0
  else:
   print("error")
 elif type==2:
  file=("/apollo/jung/bandersen/compare_full_no_parts/co2_ncep_output_"+day_time+"_24hr.txt")
  array=loadtxt(file,comments='#',delimiter=',',unpack=True)
  filepath=("/home/bandersen/iris-home/all_algorithms/image_output")
  day_label=(r+day_time+'\n')
  lat=array[0][:]
  lon=array[1][:]
  size=np.size(lat)
  if wn==1:
   co2=array[17][:]
   ncep=array[37][:]
   waver="691.875"
   coder="32"
  elif wn==2:
   co2=array[22][:]
   ncep=array[42][:]
   waver="705"
   coder="53"
  elif wn==3:
   co2=array[27][:]
   ncep=array[47][:]
   waver="715"
   coder="69"
  elif wn==4:
   co2=array[32][:]
   ncep=array[52][:]
   waver="733.125"
   coder="98"
  elif wn==5:
   co2=array[35][:]
   ncep=array[57][:]
   waver="748.125"
   coder="122"
  elif wn==6:
   co2=array[9][:]
   ncep=array[11][:]
   waver="Full"
   coder="full"
  else:
   print(error)
  lon1=[]
  lat1=[]
  lon2=[]
  lat2=[]
  lon3=[]
  lat3=[]
  lon4=[]
  lat4=[]
  if wn>=5:
   co2=np.where(co2==1,6,co2)
   co2=np.where(co2<=4,5,co2)
   for i in range(0,size):
    a=co2[i]
    b=ncep[i]
    if (a==6) and (b==0):
     c=lon[i]
     d=lat[i]
     lon1.append(c)
     lat1.append(d)
    elif (a==5) and (b==7):
     c=lon[i]
     d=lat[i]
     lon2.append(c)
     lat2.append(d)
    elif (a==6) and (b==7):
     c=lon[i]
     d=lat[i]
     lon3.append(c)
     lat3.append(d)
    elif (a==5) and (b==0):
     c=lon[i]
     d=lat[i]
     lon4.append(c)
     lat4.append(d)
    else:
     c=0
     d=0
  elif wn==6:
   co2=np.where(co2==2050,2000,co2)
   co2=np.where(co2==2100,300,co2)
   for i in range(0,size):
    a=co2[i]
    b=ncep[i]
    if (a>=900) and (b>=900):
     c=lon[i]
     d=lat[i]
     lon1.append(c)
     lat1.append(d)
    elif (a<900) and (b<900):
     c=lon[i]
     d=lat[i]
     lon2.append(c)
     lat2.append(d)
    elif (a>=900) and (b<900):
     c=lon[i]
     d=lat[i]
     lon3.append(c)
     lat3.append(d)
    elif (a<900) and (b>=900):
     c=lon[i]
     d=lat[i]
     lon4.append(c)
     lat4.append(d)
    else:
     c=0
     d=0
  else:
   print("error")
 else:
  print("error")
 fig=plt.figure()
 proj=ccrs.PlateCarree()
 ax=plt.axes(projection=proj)
 ax.add_feature(cartopy.feature.OCEAN,zorder=0,alpha=0.5,facecolor='lightskyblue')
 ax.add_feature(cartopy.feature.BORDERS,edgecolor='slategray')
 ax.add_feature(cartopy.feature.COASTLINE,edgecolor='black')
 ax.add_feature(cartopy.feature.LAND,zorder=0,edgecolor='black',alpha=0.5,facecolor='beige')
 ax.set_extent(extent,crs=proj)
 gl=ax.gridlines(crs=proj,linewidth=2,color='grey',alpha=0.5,linestyle='--',draw_labels=True)
 gl.xlocator=mticker.FixedLocator(x_ticks1)
 gl.ylocator=mticker.FixedLocator(y_ticks1)
 gl.xlabels_top=False
 gl.xlabels_bottom=False
 gl.ylabels_right=False
 gl.ylabels_left=False
 gl.xlines=True
 gl.xformatter=LONGITUDE_FORMATTER
 gl.yformatter=LATITUDE_FORMATTER
 gl.xlabel_style={'color':'red','weight':'bold','size':'40'}
 gl.ylabel_style={'color':'red','weight':'bold','size':'40'}
 g2=ax.gridlines(crs=proj,linewidth=2,color='grey',alpha=0.5,linestyle='--',draw_labels=True)
 g2.xlocator=mticker.FixedLocator(x_ticks2)
 g2.ylocator=mticker.FixedLocator(y_ticks2)
 g2.xlabels_bottom=False
 g2.ylabels_right=False
 g2.ylabels_left=True
 g2.xlines=True
 g2.xformatter=LONGITUDE_FORMATTER
 g2.yformatter=LATITUDE_FORMATTER
 g2.xlabel_style={'color':'black','weight':'bold','size':'14'}
 g2.ylabel_style={'color':'black','weight':'bold','size':'14'}
 plt.scatter(lon1,lat1,color='green',s=0.05)
 plt.title('CO2 Slicing Clear and NCEP Algorithm Clear for Wavenumber '+waver+'\n',fontsize=22,weight='bold')
 filepth=(filepath+"/confustion_co2_clr_ncep_clr_"+coder+'_'+day_time+".png")
 fig.set_size_inches(12,8.75)
 plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight')
 fig=plt.figure()
 proj=ccrs.PlateCarree()
 ax=plt.axes(projection=proj)
 ax.add_feature(cartopy.feature.OCEAN,zorder=0,alpha=0.5,facecolor='lightskyblue')
 ax.add_feature(cartopy.feature.BORDERS,edgecolor='slategray')
 ax.add_feature(cartopy.feature.COASTLINE,edgecolor='black')
 ax.add_feature(cartopy.feature.LAND,zorder=0,edgecolor='black',alpha=0.5,facecolor='beige')
 ax.set_extent(extent,crs=proj)
 gl=ax.gridlines(crs=proj,linewidth=2,color='grey',alpha=0.5,linestyle='--',draw_labels=True)
 gl.xlocator=mticker.FixedLocator(x_ticks1)
 gl.ylocator=mticker.FixedLocator(y_ticks1)
 gl.xlabels_top=False
 gl.xlabels_bottom=False
 gl.ylabels_right=False
 gl.ylabels_left=False
 gl.xlines=True
 gl.xformatter=LONGITUDE_FORMATTER
 gl.yformatter=LATITUDE_FORMATTER
 gl.xlabel_style={'color':'red','weight':'bold','size':'40'}
 gl.ylabel_style={'color':'red','weight':'bold','size':'40'}
 g2=ax.gridlines(crs=proj,linewidth=2,color='grey',alpha=0.5,linestyle='--',draw_labels=True)
 g2.xlocator=mticker.FixedLocator(x_ticks2)
 g2.ylocator=mticker.FixedLocator(y_ticks2)
 g2.xlabels_bottom=False
 g2.ylabels_right=False
 g2.ylabels_left=True
 g2.xlines=True
 g2.xformatter=LONGITUDE_FORMATTER
 g2.yformatter=LATITUDE_FORMATTER
 g2.xlabel_style={'color':'black','weight':'bold','size':'14'}
 g2.ylabel_style={'color':'black','weight':'bold','size':'14'}
 plt.scatter(lon2,lat2,color='blue',s=0.05)
 plt.title('CO2 Slicing Cloud and NCEP Algorithm Cloud for Wavenumber '+waver+'\n',fontsize=22,weight='bold')
 filepth=(filepath+"/confustion_co2_cld_ncep_cld_"+coder+'_'+day_time+".png")
 fig.set_size_inches(12,8.75)
 plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight')
 fig=plt.figure()
 proj=ccrs.PlateCarree()
 ax=plt.axes(projection=proj)
 ax.add_feature(cartopy.feature.OCEAN,zorder=0,alpha=0.5,facecolor='lightskyblue')
 ax.add_feature(cartopy.feature.BORDERS,edgecolor='slategray')
 ax.add_feature(cartopy.feature.COASTLINE,edgecolor='black')
 ax.add_feature(cartopy.feature.LAND,zorder=0,edgecolor='black',alpha=0.5,facecolor='beige')
 ax.set_extent(extent,crs=proj)
 gl=ax.gridlines(crs=proj,linewidth=2,color='grey',alpha=0.5,linestyle='--',draw_labels=True)
 gl.xlocator=mticker.FixedLocator(x_ticks1)
 gl.ylocator=mticker.FixedLocator(y_ticks1)
 gl.xlabels_top=False
 gl.xlabels_bottom=False
 gl.ylabels_right=False
 gl.ylabels_left=False
 gl.xlines=True
 gl.xformatter=LONGITUDE_FORMATTER
 gl.yformatter=LATITUDE_FORMATTER
 gl.xlabel_style={'color':'red','weight':'bold','size':'40'}
 gl.ylabel_style={'color':'red','weight':'bold','size':'40'}
 g2=ax.gridlines(crs=proj,linewidth=2,color='grey',alpha=0.5,linestyle='--',draw_labels=True)
 g2.xlocator=mticker.FixedLocator(x_ticks2)
 g2.ylocator=mticker.FixedLocator(y_ticks2)
 g2.xlabels_bottom=False
 g2.ylabels_right=False
 g2.ylabels_left=True
 g2.xlines=True
 g2.xformatter=LONGITUDE_FORMATTER
 g2.yformatter=LATITUDE_FORMATTER
 g2.xlabel_style={'color':'black','weight':'bold','size':'14'}
 g2.ylabel_style={'color':'black','weight':'bold','size':'14'}
 plt.scatter(lon3,lat3,color='red',s=0.05)
 plt.title('CO2 Slicing Clear and NCEP Algorithm Cloud for Wavenumber '+waver+'\n',fontsize=22,weight='bold')
 filepth=(filepath+"/confustion_co2_clr_ncep_cld_"+coder+'_'+day_time+".png")
 fig.set_size_inches(12,8.75)
 plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight')
 fig=plt.figure()
 proj=ccrs.PlateCarree()
 ax=plt.axes(projection=proj)
 ax.add_feature(cartopy.feature.OCEAN,zorder=0,alpha=0.5,facecolor='lightskyblue')
 ax.add_feature(cartopy.feature.BORDERS,edgecolor='slategray')
 ax.add_feature(cartopy.feature.COASTLINE,edgecolor='black')
 ax.add_feature(cartopy.feature.LAND,zorder=0,edgecolor='black',alpha=0.5,facecolor='beige')
 ax.set_extent(extent,crs=proj)
 gl=ax.gridlines(crs=proj,linewidth=2,color='grey',alpha=0.5,linestyle='--',draw_labels=True)
 gl.xlocator=mticker.FixedLocator(x_ticks1)
 gl.ylocator=mticker.FixedLocator(y_ticks1)
 gl.xlabels_top=False
 gl.xlabels_bottom=False
 gl.ylabels_right=False
 gl.ylabels_left=False
 gl.xlines=True
 gl.xformatter=LONGITUDE_FORMATTER
 gl.yformatter=LATITUDE_FORMATTER
 gl.xlabel_style={'color':'red','weight':'bold','size':'40'}
 gl.ylabel_style={'color':'red','weight':'bold','size':'40'}
 g2=ax.gridlines(crs=proj,linewidth=2,color='grey',alpha=0.5,linestyle='--',draw_labels=True)
 g2.xlocator=mticker.FixedLocator(x_ticks2)
 g2.ylocator=mticker.FixedLocator(y_ticks2)
 g2.xlabels_bottom=False
 g2.ylabels_right=False
 g2.ylabels_left=True
 g2.xlines=True
 g2.xformatter=LONGITUDE_FORMATTER
 g2.yformatter=LATITUDE_FORMATTER
 g2.xlabel_style={'color':'black','weight':'bold','size':'14'}
 g2.ylabel_style={'color':'black','weight':'bold','size':'14'}
 plt.scatter(lon4,lat4,color='purple',s=0.05)
 plt.title('CO2 Slicing Cloud and NCEP Algorithm Clear for Wavenumber '+waver+'\n',fontsize=22,weight='bold')
 filepth=(filepath+"/confustion_co2_cld_ncep_clr_"+coder+'_'+day_time+".png")
 fig.set_size_inches(12,8.75)
 plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight')

def ecmwf_ncep_count_compares(yr,mn,dy):
 b1=int(dy)
 b2=int(mn)
 b3=int(yr)
 if b1>=0 and b1<=9:
  c1="0"+str(b1)
 else:
  c1=str(b1)
 if b2>=0 and b2<=9:
  c2="0"+str(b2)
 else:
  c2=str(b2)
 c3=str(b3)
 day=c1
 mon=c2
 yer=c3
 day_time=(yer+'_'+mon+'_'+day)
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
 filepath=("/home/bandersen/iris-home/all_algorithms/image_output")
 day_label=(r+day_time+'\n')
 files_ecmwf=glob('/apollo/jung/bandersen/ecmwf_data_output_txt_'+day_time+'/*.txt')[:]
 files_ncep=glob('/apollo/jung/bandersen/ncep_data_'+day_time+'/*.nc')[:]
 sizer_e=np.size(files_ecmwf)
 sizer_n=np.size(files_ncep)
 counts_zeros_e=[]
 counts_ones_e=[]
 for i in range(0,sizer_e):
  array_a1=loadtxt(files_ecmwf[i],comments='#',delimiter=',',unpack=True)
  for ii in range(0,431):
   a1=3+ii
   zeros=np.size(np.where(array_a1[a1][:]==0))
   ones=np.size(np.where(array_a1[a1][:]==1))
   counts_zeros_e.append(zeros)
   counts_ones_e.append(ones)
 zeros_e=np.reshape(counts_zeros_e,(431,24))
 ones_e=np.reshape(counts_ones_e,(431,24))
 sum_zero_e=np.sum(zeros_e,axis=1)
 sum_ones_e=np.sum(ones_e,axis=1)
 counts_zeros_n=[]
 counts_sevens_n=[]
 for i in range(0,sizer_n):
  array_a1=cdf.Dataset(files_ncep[i])
  for ii in range(0,431):
   a1=ii
   a2=array_a1['cld_impacted'][:,a1]
   zeros=np.size(np.where(a2==0))
   sevens=np.size(np.where(a2==7))
   counts_zeros_n.append(zeros)
   counts_sevens_n.append(sevens)
 zeros_n=np.reshape(counts_zeros_n,(48,431))
 sevens_n=np.reshape(counts_sevens_n,(48,431))
 sum_zero_n=np.sum(zeros_n,axis=0)
 sum_sevens_n=np.sum(sevens_n,axis=0)
 fig,ax1=plt.subplots()
 ax1.set_xlabel('CrIS Wavenumber',color='black',fontsize=16)
 ax1.set_ylabel('FOV Count ECMWF-CAD',color='green',fontsize=16)
 ax1.scatter(wave,sum_zero_e,color='green')
 ax1.tick_params(axis='x',labelcolor='black',labelsize=16)
 ax1.tick_params(axis='y',labelcolor='green',labelsize=16)
 ax1.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1E'))
 ax2=ax1.twinx()
 ax2.set_ylabel('FOV Count NCEP Algorithm',color='red',fontsize=16)
 ax2.scatter(wave,sum_zero_n,color='red')
 ax2.tick_params(axis='y',labelcolor='red',labelsize=16)
 ax2.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1E'))
 plt.title('Count of Clear FOVs for NCEP Algorithm (red) and ECMWF-CAD (green) at all CrIS Wavenumbers'+'\n'+day_label,fontsize=22,fontweight='bold')
 filepth=(filepath+"/ecmwf_ncep_count"+day_time+".png")
 fig.set_size_inches(20,14)
 plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight')

#This produces all the observation minus analysis plots for a given day for all CO2 Slicing wavenumbers
def two_y_axis_o_a_plots(yr,mn,dy,type):
 chan_wv=[689.375,690.625,691.875,693.125,694.375,702.500,703.750,705.000,706.250,707.500,712.500,713.750,715.000,716.250,717.500,730.625,731.875,733.125,734.375,735.625,748.125,959.375]
 chan_num=[28,30,32,34,36,49,51,53,55,57,65,67,69,71,73,94,96,98,100,102,122,216]
 def counter(array):
  a=[]
  arr=-1*array
  for i in range(-300,300,1):
   b1=i*0.01
   b2=b1+0.01
   c=np.size(np.where((arr>b1)&(arr<=b2)))
   a.append(c)
  return(a)
 xs=np.arange(-3,3,0.01)
 ylab_spec=(r"FOV Count in 0.01 [K] Bin")
 def counter(array):
  a=[]
  arr=-1*array
  for i in range(-30,30,1):
   b1=i*0.1
   b2=b1+0.1
   c=np.size(np.where((arr>b1)&(arr<=b2)))
   a.append(c)
  return(a)
 xs=np.arange(-3,3,0.1)
 ylab_spec=(r"FOV Count in 0.1 [K] Bin")
 b1=int(dy)
 b2=int(mn)
 b3=int(yr)
 if b1>=0 and b1<=9:
  c1="0"+str(b1)
 else:
  c1=str(b1)
 if b2>=0 and b2<=9:
  c2="0"+str(b2)
 else:
  c2=str(b2)
 c3=str(b3)
 day=c1
 mon=c2
 yer=c3
 day_time=(yer+'_'+mon+'_'+day)
 if type==1:
  file=("/apollo/jung/bandersen/compare_full_no_parts/co2_dr_ncep_ecmwf_output_"+day_time+"_24hr.txt")
  array=loadtxt(file,comments='#',delimiter=',',unpack=True)
  filepath=("/home/bandersen/iris-home/all_algorithms/image_output")
  day_label=(r+day_time+'\n')
  forr=array[3][:]
  array=loadtxt(file,comments='#',delimiter=',',unpack=True)
  for i in range (0,22,1):
   if i<21:
    a1=17+i
    a2=38+i
    a3=60+i
    a4=82+i
    co2=array[a1][:]
    ncep=array[a2][:]
    ecmwf=array[a3][:]
    all_obs=array[a4][:]
    topo=array[2][:]
    lat=array[0][:]
    co2=np.where(co2==1,6,co2)
    co2=np.where(co2<=4,5,co2)
    co2_cld=np.asarray(counter(all_obs[np.where((co2==5))][:]))
    co2_clr=np.asarray(counter(all_obs[np.where((co2==6))][:]))
    ncep_cld=np.asarray(counter(all_obs[np.where((ncep==7))][:]))
    ncep_clr=np.asarray(counter(all_obs[np.where((ncep!=7))][:]))
    ecmwf_cld=np.asarray(counter(all_obs[np.where((ecmwf==1))][:]))
    ecmwf_clr=np.asarray(counter(all_obs[np.where((ecmwf!=1))][:]))
    fig,ax1=plt.subplots()
    ax1.set_xlabel('OBS - Analysis Brightness Temperature [K]',color='black',fontsize=10)
    ax1.set_ylabel(ylab_spec+" CO2",color='blue',fontsize=10)
    ax1.plot(xs,co2_cld,color='blue',linestyle='solid')
    ax1.tick_params(axis='y',labelcolor='blue',labelsize=10)
    ax1.tick_params(axis='x',labelcolor='black',labelsize=10)
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    ax2=ax1.twinx()
    ax2.set_ylabel(ylab_spec+" NCEP",color='red',fontsize=10)
    ax2.plot(xs,ncep_cld,color='red',linestyle='solid')
    ax2.tick_params(axis='y',labelcolor='red',labelsize=10)
    ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    plt.xlim([-3,3])
    plt.suptitle('OBS - Analysis Brightness Temperature 0.1 Degree Bins for'+'\nCO2 Slicing and NCEP Cloud Points at Wavenumber '+str(chan_wv[i]),fontsize=14,fontweight='bold')
    plt.title(day_label,fontsize=12)
    filepth=(filepath+"/oa_plot_co2_ncep_cld_"+str(chan_num[i])+"_"+day_time+".png")
    fig.set_size_inches(6.5,9)
    plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight',facecolor=(0.8392156,0.8352941,0.6627451),edgecolor='#FF5657')
    fig,ax1=plt.subplots()
    ax1.set_xlabel('OBS - Analysis Brightness Temperature [K]',color='black',fontsize=10)
    ax1.set_ylabel(ylab_spec+" CO2",color='blue',fontsize=10)
    ax1.plot(xs,co2_cld,color='blue',linestyle='solid')
    ax1.plot(xs,co2_cld,color='blue',linestyle='solid')
    ax1.tick_params(axis='x',labelcolor='black',labelsize=10)
    ax1.tick_params(axis='y',labelcolor='blue',labelsize=10)
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    ax2=ax1.twinx()
    ax2.set_ylabel(ylab_spec+" ECMWF",color='red',fontsize=10)
    ax2.plot(xs,ecmwf_cld,color='red',linestyle='solid')
    ax2.tick_params(axis='y',labelcolor='red',labelsize=10)
    ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    plt.xlim([-3,3])
    plt.suptitle('OBS - Analysis Brightness Temperature 0.1 Degree Bins for'+'\nCO2 Slicing and ECMWF Cloud Points at Wavenumber '+str(chan_wv[i]),fontsize=14,fontweight='bold')
    plt.title(day_label,fontsize=12)
    filepth=(filepath+"/oa_plot_co2_ecmwf_cld_"+str(chan_num[i])+"_"+day_time+".png")
    fig.set_size_inches(6.5,9)
    plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight',facecolor=(0.8392156,0.8352941,0.6627451),edgecolor='#FF5657')
    fig,ax1=plt.subplots()
    ax1.set_xlabel('OBS - Analysis Brightness Temperature [K]',color='black',fontsize=10)
    ax1.set_ylabel(ylab_spec+" ECMWF",color='red',fontsize=10)
    ax1.plot(xs,co2_cld,color='red',linestyle='solid')
    ax1.tick_params(axis='y',labelcolor='red',labelsize=10)
    ax1.tick_params(axis='x',labelcolor='black',labelsize=10)
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    ax2=ax1.twinx()
    ax2.set_ylabel(ylab_spec+" NCEP",color='red',fontsize=10)
    ax2.plot(xs,ncep_cld,color='red',linestyle='solid')
    ax2.tick_params(axis='y',labelcolor='red',labelsize=10)
    ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    plt.xlim([-3,3])
    plt.suptitle('OBS - Analysis Brightness Temperature 0.1 Degree Bins for'+'\nECMWF and NCEP Cloud Points at Wavenumber '+str(chan_wv[i]),fontsize=14,fontweight='bold')
    plt.title(day_label,fontsize=12)
    filepth=(filepath+"/oa_plot_ecmwf_ncep_cld_"+str(chan_num[i])+"_"+day_time+".png")
    fig.set_size_inches(6.5,9)
    plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight',facecolor=(0.8392156,0.8352941,0.6627451),edgecolor='#FF5657')
    fig,ax1=plt.subplots()
    ax1.set_xlabel('OBS - Analysis Brightness Temperature [K]',color='black',fontsize=10)
    ax1.set_ylabel(ylab_spec+" CO2",color='blue',fontsize=10)
    ax1.plot(xs,co2_clr,color='blue',linestyle='solid')
    ax1.tick_params(axis='x',labelcolor='black',labelsize=10)
    ax1.tick_params(axis='y',labelcolor='blue',labelsize=10)
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    ax2=ax1.twinx()
    ax2.set_ylabel(ylab_spec+" NCEP",color='red',fontsize=10)
    ax2.plot(xs,ncep_clr,color='red',linestyle='solid')
    ax2.tick_params(axis='y',labelcolor='red',labelsize=10)
    ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    plt.xlim([-3,3])
    plt.suptitle('OBS - Analysis Brightness Temperature 0.1 Degree Bins for'+'\nCO2 Slicing and NCEP Clear Points at Wavenumber '+str(chan_wv[i]),fontsize=14,fontweight='bold')
    plt.title(day_label,fontsize=12)
    filepth=(filepath+"/oa_plot_co2_ncep_clr_"+str(chan_num[i])+"_"+day_time+".png")
    fig.set_size_inches(6.5,9)
    plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight',facecolor=(0.8392156,0.8352941,0.6627451),edgecolor='#FF5657')
    fig,ax1=plt.subplots()
    ax1.set_xlabel('OBS - Analysis Brightness Temperature [K]',color='black',fontsize=10)
    ax1.set_ylabel(ylab_spec+" CO2",color='blue',fontsize=10)
    ax1.plot(xs,co2_clr,color='blue',linestyle='solid')
    ax1.tick_params(axis='y',labelcolor='blue',labelsize=10)
    ax1.tick_params(axis='x',labelcolor='black',labelsize=10)
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    ax2=ax1.twinx()
    ax2.set_ylabel(ylab_spec+" ECMWF",color='red',fontsize=10)
    ax2.plot(xs,ecmwf_clr,color='red',linestyle='solid')
    ax2.tick_params(axis='y',labelcolor='red',labelsize=10)
    ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    plt.xlim([-3,3])
    plt.suptitle('OBS - Analysis Brightness Temperature 0.1 Degree Bins for'+'\nCO2 Slicing and ECMWF Clear Points at Wavenumber '+str(chan_wv[i]),fontsize=14,fontweight='bold')
    plt.title(day_label,fontsize=12)
    filepth=(filepath+"/oa_plot_co2_ecmwf_clr_"+str(chan_num[i])+"_"+day_time+".png")
    fig.set_size_inches(6.5,9)
    plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight',facecolor=(0.8392156,0.8352941,0.6627451),edgecolor='#FF5657')
    fig,ax1=plt.subplots()
    ax1.set_xlabel('OBS - Analysis Brightness Temperature [K]',color='black',fontsize=10)
    ax1.set_ylabel(ylab_spec+" ECMWF",color='red',fontsize=10)
    ax1.plot(xs,co2_clr,color='red',linestyle='solid')
    ax1.tick_params(axis='y',labelcolor='red',labelsize=10)
    ax1.tick_params(axis='x',labelcolor='black',labelsize=10)
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    ax2=ax1.twinx()
    ax2.set_ylabel(ylab_spec+" NCEP",color='red',fontsize=10)
    ax2.plot(xs,ncep_clr,color='red',linestyle='solid')
    ax2.tick_params(axis='y',labelcolor='red',labelsize=10)
    ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    plt.xlim([-3,3])
    plt.suptitle('OBS - Analysis Brightness Temperature 0.1 Degree Bins for'+'\n ECMWF and NCEP Clear Points at Wavenumber '+str(chan_wv[i]),fontsize=14,fontweight='bold')
    plt.title(day_label,fontsize=12)
    filepth=(filepath+"/oa_plot_ecmwf_ncep_clr_"+str(chan_num[i])+"_"+day_time+".png")
    fig.set_size_inches(6.5,9)
    plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight',facecolor=(0.8392156,0.8352941,0.6627451),edgecolor='#FF5657')
   else:
    a1=17+i-1
    a2=38+i
    a3=60+i
    a4=82+i
    co2=array[a1][:]
    ncep=array[a2][:]
    ecmwf=array[a3][:]
    all_obs=array[a4][:]
    topo=array[2][:]
    lat=array[0][:]
    co2=np.where(co2==1,6,co2)
    co2=np.where(co2<=4,5,co2)
    co2_cld=np.asarray(counter(all_obs[np.where((co2==5))][:]))
    co2_clr=np.asarray(counter(all_obs[np.where((co2==6))][:]))
    ncep_cld=np.asarray(counter(all_obs[np.where((ncep==7))][:]))
    ncep_clr=np.asarray(counter(all_obs[np.where((ncep!=7))][:]))
    ecmwf_cld=np.asarray(counter(all_obs[np.where((ecmwf==1))][:]))
    ecmwf_clr=np.asarray(counter(all_obs[np.where((ecmwf!=1))][:]))
    fig,ax1=plt.subplots()
    ax1.set_xlabel('OBS - Analysis Brightness Temperature [K]',color='black',fontsize=10)
    ax1.set_ylabel(ylab_spec+" CO2",color='blue',fontsize=10)
    ax1.plot(xs,co2_cld,color='blue',linestyle='solid')
    ax1.tick_params(axis='y',labelcolor='blue',labelsize=10)
    ax1.tick_params(axis='x',labelcolor='black',labelsize=10)
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    ax2=ax1.twinx()
    ax2.set_ylabel(ylab_spec+" NCEP",color='red',fontsize=10)
    ax2.plot(xs,ncep_cld,color='red',linestyle='solid')
    ax2.tick_params(axis='y',labelcolor='red',labelsize=10)
    ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    plt.xlim([-3,3])
    plt.suptitle('OBS - Analysis Brightness Temperature 0.1 Degree Bins for'+'\nCO2 Slicing and NCEP Cloud Points at Wavenumber '+str(chan_wv[i]),fontsize=14,fontweight='bold')
    plt.title(day_label,fontsize=12)
    filepth=(filepath+"/oa_plot_co2_ncep_cld_"+str(chan_num[i])+"_"+day_time+".png")
    fig.set_size_inches(6.5,9)
    plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight',facecolor=(0.8392156,0.8352941,0.6627451),edgecolor='#FF5657')
    fig,ax1=plt.subplots()
    ax1.set_xlabel('OBS - Analysis Brightness Temperature [K]',color='black',fontsize=10)
    ax1.set_ylabel(ylab_spec+" CO2",color='blue',fontsize=10)
    ax1.plot(xs,co2_cld,color='blue',linestyle='solid')
    ax1.tick_params(axis='y',labelcolor='blue',labelsize=10)
    ax1.tick_params(axis='x',labelcolor='black',labelsize=10)
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    ax2=ax1.twinx()
    ax2.set_ylabel(ylab_spec+" ECMWF",color='red',fontsize=10)
    ax2.plot(xs,ecmwf_cld,color='red',linestyle='solid')
    ax2.tick_params(axis='y',labelcolor='red',labelsize=10)
    ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    plt.xlim([-3,3])
    plt.suptitle('OBS - Analysis Brightness Temperature 0.1 Degree Bins for'+'\nCO2 Slicing and ECMWF Cloud Points at Wavenumber '+str(chan_wv[i]),fontsize=14,fontweight='bold')
    plt.title(day_label,fontsize=12)
    filepth=(filepath+"/oa_plot_co2_ecmwf_cld_"+str(chan_num[i])+"_"+day_time+".png")
    fig.set_size_inches(6.5,9)
    plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight',facecolor=(0.8392156,0.8352941,0.6627451),edgecolor='#FF5657')
    fig,ax1=plt.subplots()
    ax1.set_xlabel('OBS - Analysis Brightness Temperature [K]',color='black',fontsize=10)
    ax1.set_ylabel(ylab_spec+" ECMWF",color='red',fontsize=10)
    ax1.plot(xs,co2_cld,color='red',linestyle='solid')
    ax1.tick_params(axis='y',labelcolor='red',labelsize=10)
    ax1.tick_params(axis='x',labelcolor='black',labelsize=10)
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    ax2=ax1.twinx()
    ax2.set_ylabel(ylab_spec+" NCEP",color='red',fontsize=10)
    ax2.plot(xs,ncep_cld,color='red',linestyle='solid')
    ax2.tick_params(axis='y',labelcolor='red',labelsize=10)
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    plt.xlim([-3,3])
    plt.suptitle('OBS - Analysis Brightness Temperature 0.1 Degree Bins for'+'\nECMWF and NCEP Cloud Points at Wavenumber '+str(chan_wv[i]),fontsize=14,fontweight='bold')
    plt.title(day_label,fontsize=12)
    filepth=(filepath+"/oa_plot_ecmwf_ncep_cld_"+str(chan_num[i])+"_"+day_time+".png")
    fig.set_size_inches(6.5,9)
    plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight',facecolor=(0.8392156,0.8352941,0.6627451),edgecolor='#FF5657')
    fig,ax1=plt.subplots()
    ax1.set_xlabel('OBS - Analysis Brightness Temperature [K]',color='black',fontsize=10)
    ax1.set_ylabel(ylab_spec+" CO2",color='blue',fontsize=10)
    ax1.plot(xs,co2_clr,color='blue',linestyle='solid')
    ax1.tick_params(axis='y',labelcolor='blue',labelsize=10)
    ax1.tick_params(axis='x',labelcolor='black',labelsize=10)
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    ax2=ax1.twinx()
    ax2.set_ylabel(ylab_spec+" NCEP",color='red',fontsize=10)
    ax2.plot(xs,ncep_clr,color='red',linestyle='solid')
    ax2.tick_params(axis='y',labelcolor='red',labelsize=10)
    ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    plt.xlim([-3,3])
    plt.suptitle('OBS - Analysis Brightness Temperature 0.1 Degree Bins for'+'\nCO2 Slicing and NCEP Clear Points at Wavenumber '+str(chan_wv[i]),fontsize=14,fontweight='bold')
    plt.title(day_label,fontsize=12)
    filepth=(filepath+"/oa_plot_co2_ncep_clr_"+str(chan_num[i])+"_"+day_time+".png")
    fig.set_size_inches(6.5,9)
    plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight',facecolor=(0.8392156,0.8352941,0.6627451),edgecolor='#FF5657')
    fig,ax1=plt.subplots()
    ax1.set_xlabel('OBS - Analysis Brightness Temperature [K]',color='black',fontsize=10)
    ax1.set_ylabel(ylab_spec+" CO2",color='blue',fontsize=10)
    ax1.plot(xs,co2_clr,color='blue',linestyle='solid')
    ax1.tick_params(axis='y',labelcolor='blue',labelsize=10)
    ax1.tick_params(axis='x',labelcolor='black',labelsize=10)
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    ax2=ax1.twinx()
    ax2.set_ylabel(ylab_spec+" ECMWF",color='red',fontsize=10)
    ax2.plot(xs,ecmwf_clr,color='red',linestyle='solid')
    ax2.tick_params(axis='y',labelcolor='red',labelsize=10)
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    plt.xlim([-3,3])
    plt.suptitle('OBS - Analysis Brightness Temperature 0.1 Degree Bins for'+'\nCO2 Slicing and ECMWF Clear Points at Wavenumber '+str(chan_wv[i]),fontsize=14,fontweight='bold')
    plt.title(day_label,fontsize=12)
    filepth=(filepath+"/oa_plot_co2_ecmwf_clr_"+str(chan_num[i])+"_"+day_time+".png")
    fig.set_size_inches(6.5,9)
    plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight',facecolor=(0.8392156,0.8352941,0.6627451),edgecolor='#FF5657')
    fig,ax1=plt.subplots()
    ax1.set_xlabel('OBS - Analysis Brightness Temperature [K]',color='black',fontsize=10)
    ax1.set_ylabel(ylab_spec+" ECMWF",color='red',fontsize=10)
    ax1.plot(xs,co2_clr,color='red',linestyle='solid')
    ax1.tick_params(axis='y',labelcolor='red',labelsize=10)
    ax1.tick_params(axis='x',labelcolor='black',labelsize=10)
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    ax2=ax1.twinx()
    ax2.set_ylabel(ylab_spec+" NCEP",color='red',fontsize=10)
    ax2.plot(xs,ncep_clr,color='red',linestyle='solid')
    ax2.tick_params(axis='y',labelcolor='red',labelsize=10)
    ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    plt.xlim([-3,3])
    plt.suptitle('OBS - Analysis Brightness Temperature 0.1 Degree Bins for'+'\nECMWF and NCEP Clear Points at Wavenumber '+str(chan_wv[i]),fontsize=14,fontweight='bold')
    plt.title(day_label,fontsize=12)
    filepth=(filepath+"/oa_plot_ecmwf_ncep_clr_"+str(chan_num[i])+"_"+day_time+".png")
    fig.set_size_inches(6.5,9)
    plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight',facecolor=(0.8392156,0.8352941,0.6627451),edgecolor='#FF5657')
 elif type==2:
  file=("/apollo/jung/bandersen/compare_full_no_parts/co2_ncep_output_"+day_time+"_24hr.txt")
  array=loadtxt(file,comments='#',delimiter=',',unpack=True)
  filepath=("/home/bandersen/iris-home/all_algorithms/image_output")
  day_label=(r+day_time+'\n')
  forr=array[3][:]
  array=loadtxt(file,comments='#',delimiter=',',unpack=True)
  for i in range (0,22,1):
   if i<21:
    a1=15+i
    a2=36+i
    a3=58+i
    co2=array[a1][:]
    ncep=array[a2][:]
    all_obs=array[a3][:]
    topo=array[2][:]
    lat=array[0][:]
    co2=np.where(co2==1,6,co2)
    co2=np.where(co2<=4,5,co2)
    co2_cld=np.asarray(counter(all_obs[np.where((co2==5))][:]))
    co2_clr=np.asarray(counter(all_obs[np.where((co2==6))][:]))
    ncep_cld=np.asarray(counter(all_obs[np.where((ncep==7))][:]))
    ncep_clr=np.asarray(counter(all_obs[np.where((ncep!=7))][:]))
    fig,ax1=plt.subplots()
    ax1.set_xlabel('OBS - Analysis Brightness Temperature [K]',color='black',fontsize=10)
    ax1.set_ylabel(ylab_spec+" CO2",color='blue',fontsize=10)
    ax1.plot(xs,co2_cld,color='blue',linestyle='solid')
    ax1.tick_params(axis='y',labelcolor='blue',labelsize=10)
    ax1.tick_params(axis='x',labelcolor='black',labelsize=10)
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    ax2=ax1.twinx()
    ax2.set_ylabel(ylab_spec+" NCEP",color='red',fontsize=10)
    ax2.plot(xs,ncep_cld,color='red',linestyle='solid')
    ax2.tick_params(axis='y',labelcolor='red',labelsize=10)
    ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    plt.xlim([-3,3])
    plt.suptitle('OBS - Analysis Brightness Temperature 0.1 Degree Bins for'+'\nCO2 Slicing and NCEP Cloud Points at Wavenumber '+str(chan_wv[i]),fontsize=14,fontweight='bold')
    plt.title(day_label,fontsize=12)
    filepth=(filepath+"/oa_plot_co2_ncep_cld_"+str(chan_num[i])+"_"+day_time+".png")
    fig.set_size_inches(6.5,9)
    plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight',facecolor=(0.8392156,0.8352941,0.6627451),edgecolor='#FF5657')
    fig,ax1=plt.subplots()
    ax1.set_xlabel('OBS - Analysis Brightness Temperature [K]',color='black',fontsize=10)
    ax1.set_ylabel(ylab_spec+" CO2",color='blue',fontsize=10)
    ax1.plot(xs,co2_clr,color='blue',linestyle='solid')
    ax1.tick_params(axis='x',labelcolor='black',labelsize=10)
    ax1.tick_params(axis='y',labelcolor='blue',labelsize=10)
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    ax2=ax1.twinx()
    ax2.set_ylabel(ylab_spec+" NCEP",color='red',fontsize=10)
    ax2.plot(xs,ncep_clr,color='red',linestyle='solid')
    ax2.tick_params(axis='y',labelcolor='red',labelsize=10)
    ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    plt.xlim([-3,3])
    plt.suptitle('OBS - Analysis Brightness Temperature 0.1 Degree Bins for'+'\nCO2 Slicing and NCEP Clear Points at Wavenumber '+str(chan_wv[i]),fontsize=14,fontweight='bold')
    plt.title(day_label,fontsize=12)
    filepth=(filepath+"/oa_plot_co2_ncep_clr_"+str(chan_num[i])+"_"+day_time+".png")
    fig.set_size_inches(6.5,9)
    plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight',facecolor=(0.8392156,0.8352941,0.6627451),edgecolor='#FF5657')
   else:
    a1=15+i-1
    a2=36+i
    a3=58+i
    co2=array[a1][:]
    ncep=array[a2][:]
    all_obs=array[a3][:]
    topo=array[2][:]
    lat=array[0][:]
    co2=np.where(co2==1,6,co2)
    co2=np.where(co2<=4,5,co2)
    co2_cld=np.asarray(counter(all_obs[np.where((co2==5))][:]))
    co2_clr=np.asarray(counter(all_obs[np.where((co2==6))][:]))
    ncep_cld=np.asarray(counter(all_obs[np.where((ncep==7))][:]))
    ncep_clr=np.asarray(counter(all_obs[np.where((ncep!=7))][:]))
    fig,ax1=plt.subplots()
    ax1.set_xlabel('OBS - Analysis Brightness Temperature [K]',color='black',fontsize=10)
    ax1.set_ylabel(ylab_spec+" CO2",color='blue',fontsize=10)
    ax1.plot(xs,co2_cld,color='blue',linestyle='solid')
    ax1.tick_params(axis='y',labelcolor='blue',labelsize=10)
    ax1.tick_params(axis='x',labelcolor='black',labelsize=10)
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    ax2=ax1.twinx()
    ax2.set_ylabel(ylab_spec+" NCEP",color='red',fontsize=10)
    ax2.plot(xs,ncep_cld,color='red',linestyle='solid')
    ax2.tick_params(axis='y',labelcolor='red',labelsize=10)
    ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    plt.xlim([-3,3])
    plt.suptitle('OBS - Analysis Brightness Temperature 0.1 Degree Bins for'+'\nCO2 Slicing and NCEP Cloud Points at Wavenumber '+str(chan_wv[i]),fontsize=14,fontweight='bold')
    plt.title(day_label,fontsize=12)
    filepth=(filepath+"/oa_plot_co2_ncep_cld_"+str(chan_num[i])+"_"+day_time+".png")
    fig.set_size_inches(6.5,9)
    plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight',facecolor=(0.8392156,0.8352941,0.6627451),edgecolor='#FF5657')
    fig,ax1=plt.subplots()
    ax1.set_xlabel('OBS - Analysis Brightness Temperature [K]',color='black',fontsize=10)
    ax1.set_ylabel(ylab_spec+" CO2",color='blue',fontsize=10)
    ax1.plot(xs,co2_clr,color='blue',linestyle='solid')
    ax1.tick_params(axis='y',labelcolor='blue',labelsize=10)
    ax1.tick_params(axis='x',labelcolor='black',labelsize=10)
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    ax2=ax1.twinx()
    ax2.set_ylabel(ylab_spec+" NCEP",color='red',fontsize=10)
    ax2.plot(xs,ncep_clr,color='red',linestyle='solid')
    ax2.tick_params(axis='y',labelcolor='red',labelsize=10)
    ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1E'))
    plt.xlim([-3,3])
    plt.suptitle('OBS - Analysis Brightness Temperature 0.1 Degree Bins for'+'\nCO2 Slicing and NCEP Clear Points at Wavenumber '+str(chan_wv[i]),fontsize=14,fontweight='bold')
    plt.title(day_label,fontsize=12)
    filepth=(filepath+"/oa_plot_co2_ncep_clr_"+str(chan_num[i])+"_"+day_time+".png")
    fig.set_size_inches(6.5,9)
    plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight',facecolor=(0.8392156,0.8352941,0.6627451),edgecolor='#FF5657')

def confusion_map_plots(yr,mn,dy,type,wn):
 b1=int(dy)
 b2=int(mn)
 b3=int(yr)
 if b1>=0 and b1<=9:
  c1="0"+str(b1)
 else:
  c1=str(b1)
 if b2>=0 and b2<=9:
  c2="0"+str(b2)
 else:
  c2=str(b2)
 c3=str(b3)
 day=c1
 mon=c2
 yer=c3
 day_time=(yer+'_'+mon+'_'+day)
 if type==1:
  file=("/apollo/jung/bandersen/compare_full_no_parts/co2_dr_ncep_ecmwf_output_"+day_time+"_24hr.txt")
  array=loadtxt(file,comments='#',delimiter=',',unpack=True)
  filepath=("/home/bandersen/iris-home/all_algorithms/image_output")
  day_label=(r+day_time+'\n')
  lat=array[0][:]
  topo=array[2][:]
  forr=array[3][:]
  ctpc=array[9][:]
  ctpd=array[10][:]
  ctpn=array[12][:]
  cotc=array[13][:]
  cotd=array[14][:]
  cotn=array[16][:]
  size=np.size(lat)
  ctpc=np.where(ctpc==2050,2000,ctpc)
  hig_opq_co2=np.size(np.where((ctpc>=0)&(ctpc<440)&(cotc>0.95))[:])
  hig_thk_co2=np.size(np.where((ctpc>=0)&(ctpc<440)&(cotc>=0.5)&(cotc<=0.95))[:])
  hig_thn_co2=np.size(np.where((ctpc>=0)&(ctpc<440)&(cotc<0.5))[:])
  mid_opq_co2=np.size(np.where((ctpc>=440)&(ctpc<660)&(cotc>0.95))[:])
  mid_thk_co2=np.size(np.where((ctpc>=440)&(ctpc<660)&(cotc>=0.5)&(cotc<=0.95))[:])
  mid_thn_co2=np.size(np.where((ctpc>=440)&(ctpc<660)&(cotc<0.5))[:])
  low_co2=np.size(np.where((ctpc>=660)&(ctpc<900))[:])
  clr_co2=np.size(np.where(ctpc==2000)[:])
  sfc_co2=np.size(np.where(ctpc==2100)[:])
  error_co2=np.size(np.where(ctpc==-7777)[:])
  hig_opq_nc=np.size(np.where((ctpn>=0)&(ctpn<440)&(cotn>0.95)&(cotn<1.01))[:])
  hig_thk_nc=np.size(np.where((ctpn>=0)&(ctpn<440)&(cotn>=0.5)&(cotn<=0.95))[:])
  hig_thn_nc=np.size(np.where((ctpn>=0)&(ctpn<440)&(cotn>=0)&(cotn<0.5))[:])
  mid_opq_nc=np.size(np.where((ctpn>=440)&(ctpn<660)&(cotn>0.95)&(cotn<1.01))[:])
  mid_thk_nc=np.size(np.where((ctpn>=440)&(ctpn<660)&(cotn>=0.5)&(cotn<=0.95))[:])
  mid_thn_nc=np.size(np.where((ctpn>=440)&(ctpn<660)&(cotn>=0)&(cotn<0.5))[:])
  low_nc=np.size(np.where((ctpn>=660)&(ctpn<900))[:])
  sfc_nc=np.size(np.where((ctpn>=900)&(ctpn<1010))[:])
  error_nc=np.size(np.where((ctpn<-2000))[:])
  clr_nc=int(size-hig_opq_nc-mid_opq_nc-hig_thk_nc-mid_thk_nc-hig_thn_nc-mid_thn_nc-low_nc-error_nc-sfc_nc)
  hig_opq_dr=np.size(np.where((ctpd>=0)&(ctpd<440)&(cotd<=5)&(cotd>3.75))[:])
  hig_thk_dr=np.size(np.where((ctpd>=0)&(ctpd<440)&(cotd<=3.75)&(cotd>=1))[:])
  hig_thn_dr=np.size(np.where((ctpd>=0)&(ctpd<440)&(cotd<=1)&(cotd>=0))[:])
  mid_opq_dr=np.size(np.where((ctpd>=440)&(ctpd<660)&(cotd<=5)&(cotd>3.75))[:])
  mid_thk_dr=np.size(np.where((ctpd>=440)&(ctpd<660)&(cotd<=3.75)&(cotd>=1))[:])
  mid_thn_dr=np.size(np.where((ctpd>=440)&(ctpd<660)&(cotd<=1)&(cotd>=0))[:])
  low_dr=np.size(np.where((ctpd>=660)&(ctpd<900))[:])
  sfc_dr=np.size(np.where((ctpd>=900)&(ctpd<1010))[:])
  error_dr=np.size(np.where((ctpd<-2000))[:])
  clr_dr=int(size-hig_opq_dr-mid_opq_dr-hig_thk_dr-mid_thk_dr-hig_thn_dr-mid_thn_dr-low_dr-error_dr)
  labels_co2=['High (<440 hPa)'+'\n'+r'Opaque ($ N \epsilon > 0.95$)','High (<440 hPa)'+'\n'+r'Thick ($0.5 \geq N \epsilon \geq 0.95$)','High (<440 hPa)'+'\n'+r'Thin ($ N \epsilon < 0.5$)','Middle (660-44$
  co2=np.array([hig_opq_co2,hig_thk_co2,hig_thn_co2,mid_opq_co2,mid_thk_co2,mid_thn_co2,low_co2,clr_co2,sfc_co2,error_co2])
  bcolors=['red','pink','orange','greenyellow','seagreen','darkgreen','blue','purple','black','grey']
  fig,ax=plt.subplots()
  ax.pie(co2,labels=labels_co2,colors=bcolors,autopct='%1.1f%%',shadow=True,startangle=90,textprops={'fontsize':14})
  ax.axis('equal')
  plt.title('Cloud Category for CO2 Slicing '+'\n'+day_label,fontsize=22,fontweight='bold')
  filepth=(filepath+"/pie_chart_co2_"+day_time+".png")
  fig.set_size_inches(10,10)
  plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight')
  ncep=np.array([hig_opq_nc,hig_thk_nc,hig_thn_nc,mid_opq_nc,mid_thk_nc,mid_thn_nc,low_nc,clr_nc,sfc_nc,error_nc])
  bcolors=['red','pink','orange','greenyellow','seagreen','darkgreen','blue','purple','black','grey']
  fig,ax=plt.subplots()
  ax.pie(ncep,labels=labels_co2,colors=bcolors,autopct='%1.1f%%',shadow=True,startangle=90,textprops={'fontsize':14})
  ax.axis('equal')
  plt.title('Cloud Category for NCEP Algorithm '+'\n'+day_label,fontsize=22,fontweight='bold')
  filepth=(filepath+"/pie_chart_ncep_"+day_time+".png")
  fig.set_size_inches(10,10)
  plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight')
  dr=np.array([hig_opq_dr,hig_thk_dr,hig_thn_dr,mid_opq_dr,mid_thk_dr,mid_thn_dr,low_dr,clr_dr,sfc_dr,error_dr])
  fig,ax=plt.subplots()
  ax.pie(dr,labels=labels_co2,colors=bcolors,autopct='%1.1f%%',shadow=True,startangle=90,textprops={'fontsize':14})
  ax.axis('equal')
  plt.title('Cloud Category for Dual Regression '+'\n'+day_label,fontsize=22,fontweight='bold')
  filepth=(filepath+"/pie_chart_dr_"+day_time+".png")
  fig.set_size_inches(10,10)
  plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight')
 elif type==2:
  file=("/apollo/jung/bandersen/compare_full_no_parts/co2_ncep_output_"+day_time+"_24hr.txt")
  array=loadtxt(file,comments='#',delimiter=',',unpack=True)
  filepath=("/home/bandersen/iris-home/all_algorithms/image_output")
  day_label=(r+day_time+'\n')
  forr=array[3][:]
  array=loadtxt(file,comments='#',delimiter=',',unpack=True)
  lat=array[0][:]
  topo=array[2][:]
  forr=array[3][:]
  ctpc=array[9][:]
  ctpn=array[11][:]
  cotc=array[12][:]
  cotn=array[14][:]
  size=np.size(lat)
  ctpc=np.where(ctpc==2050,2000,ctpc)
  hig_opq_co2=np.size(np.where((ctpc>=0)&(ctpc<440)&(cotc>0.95))[:])
  hig_thk_co2=np.size(np.where((ctpc>=0)&(ctpc<440)&(cotc>=0.5)&(cotc<=0.95))[:])
  hig_thn_co2=np.size(np.where((ctpc>=0)&(ctpc<440)&(cotc<0.5))[:])
  mid_opq_co2=np.size(np.where((ctpc>=440)&(ctpc<660)&(cotc>0.95))[:])
  mid_thk_co2=np.size(np.where((ctpc>=440)&(ctpc<660)&(cotc>=0.5)&(cotc<=0.95))[:])
  mid_thn_co2=np.size(np.where((ctpc>=440)&(ctpc<660)&(cotc<0.5))[:])
  low_co2=np.size(np.where((ctpc>=660)&(ctpc<900))[:])
  clr_co2=np.size(np.where(ctpc==2000)[:])
  sfc_co2=np.size(np.where(ctpc==2100)[:])
  error_co2=np.size(np.where(ctpc==-7777)[:])
  hig_opq_nc=np.size(np.where((ctpn>=0)&(ctpn<440)&(cotn>0.95)&(cotn<1.01))[:])
  hig_thk_nc=np.size(np.where((ctpn>=0)&(ctpn<440)&(cotn>=0.5)&(cotn<=0.95))[:])
  hig_thn_nc=np.size(np.where((ctpn>=0)&(ctpn<440)&(cotn>=0)&(cotn<0.5))[:])
  mid_opq_nc=np.size(np.where((ctpn>=440)&(ctpn<660)&(cotn>0.95)&(cotn<1.01))[:])
  mid_thk_nc=np.size(np.where((ctpn>=440)&(ctpn<660)&(cotn>=0.5)&(cotn<=0.95))[:])
  mid_thn_nc=np.size(np.where((ctpn>=440)&(ctpn<660)&(cotn>=0)&(cotn<0.5))[:])
  low_nc=np.size(np.where((ctpn>=660)&(ctpn<900))[:])
  sfc_nc=np.size(np.where((ctpn>=900)&(ctpn<1010))[:])
  error_nc=np.size(np.where((ctpn<-2000))[:])
  clr_nc=int(size-hig_opq_nc-mid_opq_nc-hig_thk_nc-mid_thk_nc-hig_thn_nc-mid_thn_nc-low_nc-error_nc-sfc_nc)
  labels_co2=['High (<440 hPa)'+'\n'+r'Opaque ($ N \epsilon > 0.95$)','High (<440 hPa)'+'\n'+r'Thick ($0.5 \geq N \epsilon \geq 0.95$)','High (<440 hPa)'+'\n'+r'Thin ($ N \epsilon < 0.5$)','Middle (660-44$
  co2=np.array([hig_opq_co2,hig_thk_co2,hig_thn_co2,mid_opq_co2,mid_thk_co2,mid_thn_co2,low_co2,clr_co2,sfc_co2,error_co2])
  bcolors=['red','pink','orange','greenyellow','seagreen','darkgreen','blue','purple','black','grey']
  fig,ax=plt.subplots()
  ax.pie(co2,labels=labels_co2,colors=bcolors,autopct='%1.1f%%',shadow=True,startangle=90,textprops={'fontsize':14})
  ax.axis('equal')
  plt.title('Cloud Category for CO2 Slicing '+'\n'+day_label,fontsize=22,fontweight='bold')
  filepth=(filepath+"/pie_chart_co2_"+day_time+".png")
  fig.set_size_inches(10,10)
  plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight')
  ncep=np.array([hig_opq_nc,hig_thk_nc,hig_thn_nc,mid_opq_nc,mid_thk_nc,mid_thn_nc,low_nc,clr_nc,sfc_nc,error_nc])
  bcolors=['red','pink','orange','greenyellow','seagreen','darkgreen','blue','purple','black','grey']
  fig,ax=plt.subplots()
  ax.pie(ncep,labels=labels_co2,colors=bcolors,autopct='%1.1f%%',shadow=True,startangle=90,textprops={'fontsize':14})
  ax.axis('equal')
  plt.title('Cloud Category for NCEP Algorithm '+'\n'+day_label,fontsize=22,fontweight='bold')
  filepth=(filepath+"/pie_chart_ncep_"+day_time+".png")
  fig.set_size_inches(10,10)
  plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight')
