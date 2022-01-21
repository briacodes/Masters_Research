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
filepath=("/home/bandersen/iris-home")
extent=[-180,180,-90,90]
x_ticks1=[-180,-165,-150,-135,-120,-105,-90,-75,-60,-45,-30,-15,0,15,30,45,60,75,90,105,120,135,150,165,180]
y_ticks1=[-90,-75,-60,-45,-30,-15,0,15,30,45,60,75,90]
x_ticks2=[-180,-135,-90,-45,0,45,90,135,180]
y_ticks2=[-90,-60,-30,0,30,60,90]
cmap2=colors.ListedColormap(['red','orange','green','blue','white'])
boundaries2=[0,1,2,3,4,5,6]
norm2=colors.BoundaryNorm(boundaries2,cmap2.N,clip=True)
def lev_cat_sel(i,array):
 a=array[i]
 if a>=0 and a<32:
  b=1
 elif a>=32 and a<63:
  b=2
 elif a>=63 and a<94:
  b=3
 elif a>=94 and a<125:
  b=4
 elif a>=125:
  b=5
 else:
  b=6
 return(b)

h_a=0.5
ch_a="5mb"
hh_a=str('%.1f'%h_a)
day_time_a=("2020_08_21_03_18_"+hh_a+"_"+ch_a)
day_label_a=(r'$\bf{Date}$ 2020-08-21 $\bf{Time}$ 03:00-18:00 $\bf{A(\nu)\geq}$ '+hh_a+r' $\bf{Chan}$ 4X5'+r' Minus Bias'+'\n'+'\n'+'\n')
file_a=("/apollo/jung/bandersen/compare_full_no_parts/compare_output_2020_08_21_03_18_chan_"+ch_a+"_h_"+hh_a+".txt")
array_a=loadtxt(file_a,comments='#',delimiter=',',unpack=True)
lat_a=array_a[0][:]
lon_a=array_a[1][:]
lev_a=array_a[9][:]
size_a=np.size(lat_a)
lev_cat_a=[]
for a in range(0,size_a):
 b=lev_cat_sel(a,lev_a)
 lev_cat_a.append(b)
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
g2.xlabel_style={'color':'white','weight':'bold','size':'16'}
g2.ylabel_style={'color':'white','weight':'bold','size':'16'}
data=plt.scatter(lon_a,lat_a,c=lev_cat_a,cmap=cmap2,s=2)
plt.suptitle('      CO2 Slicing Channesl Selected',weight='bold',fontsize=35)
plt.title(day_label_a,fontsize=18)
cb=plt.colorbar(data,aspect=30,fraction=.045,pad=.04)
cb.set_ticks([1.2,2.0,2.8,3.6,4.4])
cb.set_ticklabels(['691.875 Cluster','705.000 Cluster','715.000 Cluster','733.125 Cluster','CLR'])
cb.ax.tick_params(labelsize=25)
filepth=(filepath+"/LEV_CO2_cat_"+day_time_a+".png")
fig.set_size_inches(11.75,9.25)
plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight',facecolor='#AEC0C9',edgecolor='#57575F')

h_b=1.0
ch_b="5mb"
hh_b=str('%.1f'%h_b)
day_time_b=("2020_08_21_03_18_"+hh_b+"_"+ch_b)
day_label_b=(r'$\bf{Date}$ 2020-08-21 $\bf{Time}$ 03:00-18:00 $\bf{A(\nu)\geq}$ '+hh_b+r' $\bf{Chan}$ 4X5'+r' Minus Bias'+'\n'+'\n'+'\n')
file_b=("/apollo/jung/bandersen/compare_full_no_parts/compare_output_2020_08_21_03_18_chan_"+ch_b+"_h_"+hh_b+".txt")
array_b=loadtxt(file_b,comments='#',delimiter=',',unpack=True)
lat_b=array_b[0][:]
lon_b=array_b[1][:]
lev_b=array_b[9][:]
size_b=np.size(lat_b)
lev_cat_b=[]
for a in range(0,size_b):
 b=lev_cat_sel(a,lev_b)
 lev_cat_b.append(b)
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
g2.xlabel_style={'color':'white','weight':'bold','size':'16'}
g2.ylabel_style={'color':'white','weight':'bold','size':'16'}
data=plt.scatter(lon_b,lat_b,c=lev_cat_b,cmap=cmap2,s=2)
plt.suptitle('      CO2 Slicing Channesl Selected',weight='bold',fontsize=35)
plt.title(day_label_b,fontsize=18)
cb=plt.colorbar(data,aspect=30,fraction=.045,pad=.04)
cb.set_ticks([1.2,2.0,2.8,3.6,4.4])
cb.set_ticklabels(['691.875 Cluster','705.000 Cluster','715.000 Cluster','733.125 Cluster','CLR'])
cb.ax.tick_params(labelsize=25)
filepth=(filepath+"/LEV_CO2_cat_"+day_time_b+".png")
fig.set_size_inches(11.75,9.25)
plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight',facecolor='#AEC0C9',edgecolor='#57575F')

h_z=0.5
ch_z="5nb"
hh_z=str('%.1f'%h_z)
day_time_z=("2020_08_21_03_18_"+hh_z+"_"+ch_z)
day_label_z=(r'$\bf{Date}$ 2020-08-21 $\bf{Time}$ 03:00-18:00 $\bf{A(\nu)\geq}$ '+hh_z+r' $\bf{Chan}$ 4X5'+r' No Bias'+'\n'+'\n'+'\n')
file_z=("/apollo/jung/bandersen/compare_full_no_parts/compare_output_2020_08_21_03_18_chan_"+ch_z+"_h_"+hh_z+".txt")
array_z=loadtxt(file_z,comments='#',delimiter=',',unpack=True)
lat_z=array_z[0][:]
lon_z=array_z[1][:]
lev_z=array_z[9][:]
size_z=np.size(lat_z)
lev_cat_z=[]
for a in range(0,size_z):
 b=lev_cat_sel(a,lev_z)
 lev_cat_z.append(b)
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
gl.xlabels_zottom=False
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
g2.xlabels_zottom=False
g2.ylabels_right=False
g2.ylabels_left=True
g2.xlines=True
g2.xformatter=LONGITUDE_FORMATTER
g2.yformatter=LATITUDE_FORMATTER
g2.xlabel_style={'color':'white','weight':'bold','size':'16'}
g2.ylabel_style={'color':'white','weight':'bold','size':'16'}
data=plt.scatter(lon_z,lat_z,c=lev_cat_z,cmap=cmap2,s=2)
plt.suptitle('      CO2 Slicing Channesl Selected',weight='bold',fontsize=35)
plt.title(day_label_z,fontsize=18)
cb=plt.colorbar(data,aspect=30,fraction=.045,pad=.04)
cb.set_ticks([1.2,2.0,2.8,3.6,4.4])
cb.set_ticklabels(['691.875 Cluster','705.000 Cluster','715.000 Cluster','733.125 Cluster','CLR'])
cb.ax.tick_params(labelsize=25)
filepth=(filepath+"/LEV_CO2_cat_"+day_time_z+".png")
fig.set_size_inches(11.75,9.25)
plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight',facecolor='#AEC0C9',edgecolor='#57575F')

h_y=1.0
ch_y="5nb"
hh_y=str('%.1f'%h_y)
day_time_y=("2020_08_21_03_18_"+hh_y+"_"+ch_y)
day_label_y=(r'$\bf{Date}$ 2020-08-21 $\bf{Time}$ 03:00-18:00 $\bf{A(\nu)\geq}$ '+hh_y+r' $\bf{Chan}$ 4X5'+r' No Bias'+'\n'+'\n'+'\n')
file_y=("/apollo/jung/bandersen/compare_full_no_parts/compare_output_2020_08_21_03_18_chan_"+ch_y+"_h_"+hh_y+".txt")
array_y=loadtxt(file_y,comments='#',delimiter=',',unpack=True)
lat_y=array_y[0][:]
lon_y=array_y[1][:]
lev_y=array_y[9][:]
size_y=np.size(lat_y)
lev_cat_y=[]
for a in range(0,size_y):
 b=lev_cat_sel(a,lev_y)
 lev_cat_y.append(b)
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
gl.xlabels_yottom=False
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
g2.xlabels_yottom=False
g2.ylabels_right=False
g2.ylabels_left=True
g2.xlines=True
g2.xformatter=LONGITUDE_FORMATTER
g2.yformatter=LATITUDE_FORMATTER
g2.xlabel_style={'color':'white','weight':'bold','size':'16'}
g2.ylabel_style={'color':'white','weight':'bold','size':'16'}
data=plt.scatter(lon_y,lat_y,c=lev_cat_y,cmap=cmap2,s=2)
plt.suptitle('      CO2 Slicing Channesl Selected',weight='bold',fontsize=35)
plt.title(day_label_y,fontsize=18)
cb=plt.colorbar(data,aspect=30,fraction=.045,pad=.04)
cb.set_ticks([1.2,2.0,2.8,3.6,4.4])
cb.set_ticklabels(['691.875 Cluster','705.000 Cluster','715.000 Cluster','733.125 Cluster','CLR'])
cb.ax.tick_params(labelsize=25)
filepth=(filepath+"/LEV_CO2_cat_"+day_time_y+".png")
fig.set_size_inches(11.75,9.25)
plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight',facecolor='#AEC0C9',edgecolor='#57575F')

h_x=0.5
ch_x="5pb"
hh_x=str('%.1f'%h_x)
day_time_x=("2020_08_21_03_18_"+hh_x+"_"+ch_x)
day_label_x=(r'$\bf{Date}$ 2020-08-21 $\bf{Time}$ 03:00-18:00 $\bf{A(\nu)\geq}$ '+hh_x+r' $\bf{Chan}$ 4X5'+r' Plus Bias'+'\n'+'\n'+'\n')
file_x=("/apollo/jung/bandersen/compare_full_no_parts/compare_output_2020_08_21_03_18_chan_"+ch_x+"_h_"+hh_x+".txt")
array_x=loadtxt(file_x,comments='#',delimiter=',',unpack=True)
lat_x=array_x[0][:]
lon_x=array_x[1][:]
lev_x=array_x[9][:]
size_x=np.size(lat_x)
lev_cat_x=[]
for a in range(0,size_x):
 b=lev_cat_sel(a,lev_x)
 lev_cat_x.append(b)
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
gl.xlabels_xottom=False
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
g2.xlabels_xottom=False
g2.ylabels_right=False
g2.ylabels_left=True
g2.xlines=True
g2.xformatter=LONGITUDE_FORMATTER
g2.yformatter=LATITUDE_FORMATTER
g2.xlabel_style={'color':'white','weight':'bold','size':'16'}
g2.ylabel_style={'color':'white','weight':'bold','size':'16'}
data=plt.scatter(lon_x,lat_x,c=lev_cat_x,cmap=cmap2,s=2)
plt.suptitle('      CO2 Slicing Channesl Selected',weight='bold',fontsize=35)
plt.title(day_label_x,fontsize=18)
cb=plt.colorbar(data,aspect=30,fraction=.045,pad=.04)
cb.set_ticks([1.2,2.0,2.8,3.6,4.4])
cb.set_ticklabels(['691.875 Cluster','705.000 Cluster','715.000 Cluster','733.125 Cluster','CLR'])
cb.ax.tick_params(labelsize=25)
filepth=(filepath+"/LEV_CO2_cat_"+day_time_x+".png")
fig.set_size_inches(11.75,9.25)
plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight',facecolor='#AEC0C9',edgecolor='#57575F')

h_w=1.0
ch_w="5pb"
hh_w=str('%.1f'%h_w)
day_time_w=("2020_08_21_03_18_"+hh_w+"_"+ch_w)
day_label_w=(r'$\bf{Date}$ 2020-08-21 $\bf{Time}$ 03:00-18:00 $\bf{A(\nu)\geq}$ '+hh_w+r' $\bf{Chan}$ 4X5'+r' Plus Bias'+'\n'+'\n'+'\n')
file_w=("/apollo/jung/bandersen/compare_full_no_parts/compare_output_2020_08_21_03_18_chan_"+ch_w+"_h_"+hh_w+".txt")
array_w=loadtxt(file_w,comments='#',delimiter=',',unpack=True)
lat_w=array_w[0][:]
lon_w=array_w[1][:]
lev_w=array_w[9][:]
size_w=np.size(lat_w)
lev_cat_w=[]
for a in range(0,size_w):
 b=lev_cat_sel(a,lev_w)
 lev_cat_w.append(b)
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
gl.xlabels_wottom=False
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
g2.xlabels_wottom=False
g2.ylabels_right=False
g2.ylabels_left=True
g2.xlines=True
g2.xformatter=LONGITUDE_FORMATTER
g2.yformatter=LATITUDE_FORMATTER
g2.xlabel_style={'color':'white','weight':'bold','size':'16'}
g2.ylabel_style={'color':'white','weight':'bold','size':'16'}
data=plt.scatter(lon_w,lat_w,c=lev_cat_w,cmap=cmap2,s=2)
plt.suptitle('      CO2 Slicing Channesl Selected',weight='bold',fontsize=35)
plt.title(day_label_w,fontsize=18)
cb=plt.colorbar(data,aspect=30,fraction=.045,pad=.04)
cb.set_ticks([1.2,2.0,2.8,3.6,4.4])
cb.set_ticklabels(['691.875 Cluster','705.000 Cluster','715.000 Cluster','733.125 Cluster','CLR'])
cb.ax.tick_params(labelsize=25)
filepth=(filepath+"/LEV_CO2_cat_"+day_time_w+".png")
fig.set_size_inches(11.75,9.25)
plt.savefig(filepth,dpi=fig.dpi,bbox_inches='tight',facecolor='#AEC0C9',edgecolor='#57575F')

print("Done")
