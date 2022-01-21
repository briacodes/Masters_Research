#Run before in terminal: module load license_intel intel hdf hdf5 netcdf4 miniconda

###   This code runs the full CO2 slicing technique then outputs a txt file with the information.  	###
###   It also provides the colocated VIIRS data in the output file.                 			###

import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt
import numpy as np
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
from scipy import stats
import netCDF4 as cdf

def find_nearest(array,value):
 array=np.asarray(array)
 idx=(np.abs(array-value)).argmin()
 return array[idx]

def find_nearest_index(array,value):
 array=np.asarray(array)
 idx=(np.abs(array-value)).argmin()
 return(idx)

def pbl(num,pr,te):
 p=pr[num]
 t=te[num]
 a=p[np.where(np.diff(t)>=-0.5)[0]+1]
 try:
  b=int(np.max(a[np.where(a>=850)]))
  c=find_nearest_index(p,b)
  return(c)
 except ValueError:
  return(find_nearest_index(p,850))

def trop(num,pr,te,da,la):
 p=pr[num][:]
 t=te[num][:]
 d=[int(i)  for  i  in  str(da[num])]
 if (d[4]==1 and d[5]==1 and d[6]>=1) or (d[4]==1 and d[5]==2) or (d[4]==0 and d[5]==1) or (d[4]==0 and d[5]==2) or (d[4]==0 and d[5]==3 and d[6]<=1):
  e=1  #NH  Winter  as  Mid  November  to  Mid  March
 elif (d[4]==0 and d[5]==5 and d[6]>=1) or (d[4]==0 and d[5]==6) or (d[4]==0 and d[5]==7) or (d[4]==0 and d[5]==8) or (d[4]==0 and d[5]==9 and d[6]<=1):
  e=2  #NH  Summer  as  Mid  May  to  Mid  September
 else:
  e=3  #Spring and Fall  Mid  March  to  Mid  May and Mid  September  to  Mid  November
 l=la[num]
 if l>=60 and e==1:
  a=p[np.where(np.diff(t)>=-1)[0]+1]
  b=(np.max(a[np.where(a<=325)]))
  c=find_nearest_index(p,b)
 elif l<=-60 and e==1:
  a=p[np.where(np.diff(t)>=-1)[0]+1]
  b=(np.max(a[np.where(a<=350)]))
  c=find_nearest_index(p,b)
 elif l<=-60 and e==2:
  a=p[np.where(np.diff(t)>=-1)[0]+1]
  b=(np.max(a[np.where(a<=375)]))
  c=find_nearest_index(p,b)
 elif l>=60 and e==2:
  a=p[np.where(np.diff(t)>=-1)[0]+1]
  b=(np.max(a[np.where(a<=325)]))
  c=find_nearest_index(p,b)
 elif l>=60 and e==3:
  a=p[np.where(np.diff(t)>=-1)[0]+1]
  b=(np.max(a[np.where(a<=375)]))
  c=find_nearest_index(p,b)
 elif l<=-60 and e==3:
  a=p[np.where(np.diff(t)>=-1)[0]+1]
  b=(np.max(a[np.where(a<=325)]))
  c=find_nearest_index(p,b)
 else:
  a=p[np.where(np.diff(t)>=-1)[0]+1]
  b=(np.max(a[np.where(a<=250)]))
  c=find_nearest_index(p,b)
 return(c)

def full_thing(st1,st2,sp1,sp2,dt1,dt2,mon1,ye1,hs):
 h1=hs
 def h2p(height):
  if height<=44307:
   return(round((1013.25)*((1-((2.25694436*(10**(-5)))*(height)))**(5.2553026)),1))
  else:
   return(0)
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

 def planck_wn(wn,bt):
  c_1=1.191042*(10**-5)
  c_2=1.4387752
  a=np.exp(c_2*wn/bt)
  b=wn**3
  c=c_1*b/(a-1)
  return(c)

 def find_nearest(array,value):
  array=np.asarray(array)
  idx=(np.abs(array-value)).argmin()
  return array[idx]

 def find_nearest_index(array,value):
  array=np.asarray(array)
  idx=(np.abs(array-value)).argmin()
  return(idx)

 def integrate_one(a,wave1,ob,cr,bi):
  d=wave[wave1]
  bbb=bi[wave1]
  b1=ob[a,wave1]
  c1=cr[a,wave1]+bbb
  b2=planck_wn(d,b1)
  c2=planck_wn(d,c1)
  e=c2-b2
  return(e)

 def integrate_two(a,wave1,l1,l2,ta,te):
  e=wave[wave1]
  b1=ta[a,wave1,l1]
  c1=te[a,l1]
  d1=planck_wn(e,c1)
  b2=ta[a,wave1,l2]
  c2=te[a,l2]
  d2=planck_wn(e,c2)
  g=b2*d2-b1*d1
  return(g)

 def integrate_five(a,wave1,wave2,ob,cr,ta,te,pr,da,la,bi):  #where  wave1  is  small
  pc=trop(a,pr,te,da,la)-1
  ps=pbl(a,pr,te)+1
  arrayc1=[]
  arrayc2=[]
  b1=integrate_one(a,wave1,ob,cr,bi)
  b2=integrate_one(a,wave2,ob,cr,bi)
  b=b1/b2
# for i in range(pc,ps,1):
  for i in range(pc,ps,-1):
   c1=integrate_two(a,wave1,i,i-1,ta,te)
   c2=integrate_two(a,wave2,i,i-1,ta,te)
   arrayc1.append(c1)
   arrayc2.append(c2)
  sum1=np.cumsum(arrayc1[::-1])
  sum2=np.cumsum(arrayc2[::-1])
  size=np.size(sum1)
  c1=sum1[::-1]
  c2=sum2[::-1]
  results=[]
  for ii in range(0,size):
   c=(c1[ii]/c2[ii])-b
   results.append(c)
  return(results)

 def integrate_six(a,pr,te,da,la):
  pc=trop(a,pr,te,da,la)+1
  ps=pbl(a,pr,te)-1
  array=[]
  for i in range(pc,ps,-1):
   c=pr[a,i]
   array.append(c)
  b=array[::-1]
  d=b[::-1]
  return(d)

 def co2(a,num,cr,ob,bi):
  b=cr[a,num]
  c=ob[a,num]
  d=bi[num]
  e=b-c+d
  return(e)

 def trim_ave(a):
  b=np.mean(a)
  c=np.std(a)
  d=abs(a-b)
  e=np.where(d<c)[0][:]
  g=np.size(e)
  h=[]
  for i in range(0,g):
   j=e[i]
   k=a[j]
   h.append(k)
  l=np.mean(h)
  return(l)

 def optics(a,ob,cr,te,bi,pr,pres):
  b=find_nearest_index(wave,960)
  c=find_nearest_index(pr[a],pres)
  d=ob[a,b]
  e=cr[a,b]
  g=te[a,c]
  h=bi[b]
  j=(d-e-h)/(g-e-h)
  return(j)

 def optics_1(a,ob,cr,te,bi):
  b=find_nearest_index(wave,960)
  c=ob[a,b]
  d=cr[a,b]
  e=te[a,0]
  g=bi[b]
  h=(c-d-g)/(e-d-g)
  return(h)

 def ir_test_ctp(a,ob,te,pr):
  b=find_nearest_index(wave,960)
  c=ob[a,b]
  d=te[a][:]
  e=find_nearest_index(d,c)
  g=te[a][e]
  h=np.abs(c-g)
  if h<=0.5:
   i=pr[a][e]
  else:
   i=-7777
  return(i)

 def level(a,cr,ob,bi,te):
  hh=h1
  b1=co2(a,28,cr,ob,bi)	#689.375
  b2=co2(a,30,cr,ob,bi)	#690.625
  b3=co2(a,32,cr,ob,bi)	#691.875
  b4=co2(a,34,cr,ob,bi)	#693.125
  b5=co2(a,36,cr,ob,bi)	#694.375
  if b1>=hh and b2>=hh and b3>=hh and b4>=hh and b5>=hh:
   g=1
  elif b1>=hh and b2>=hh and b3>=hh and b4>=hh and b5<hh:
   g=2
  elif b1>=hh and b2>=hh and b3>=hh and b4<hh and b5>=hh:
   g=3
  elif b1>=hh and b2>=hh and b3<hh and b4>=hh and b5>=hh:
   g=4
  elif b1>=hh and b2<hh and b3>=hh and b4>=hh and b5>=hh:
   g=5
  elif b1<hh and b2>=hh and b3>=hh and b4>=hh and b5>=hh:
   g=6
  elif b1>=hh and b2>=hh and b3>=hh and b4<hh and b5<hh:
   g=7
  elif b1>=hh and b2>=hh and b3<hh and b4>=hh and b5<hh:
   g=8
  elif b1>=hh and b2<hh and b3>=hh and b4>=hh and b5<hh:
   g=9
  elif b1<hh and b2>=hh and b3>=hh and b4>=hh and b5<hh:
   g=10
  elif b1>=hh and b2>=hh and b3<hh and b4<hh and b5>=hh:
   g=11
  elif b1>=hh and b2<hh and b3>=hh and b4<hh and b5>=hh:
   g=12
  elif b1<hh and b2>=hh and b3>=hh and b4<hh and b5>=hh:
   g=13
  elif b1>=hh and b2<hh and b3<hh and b4>=hh and b5>=hh:
   g=14
  elif b1<hh and b2>=hh and b3<hh and b4>=hh and b5>=hh:
   g=15
  elif b1<hh and b2<hh and b3>=hh and b4>=hh and b5>=hh:
   g=16
  elif b1>=hh and b2>=hh and b3<hh and b4<hh and b5<hh:
   g=17
  elif b1>=hh and b2<hh and b3>=hh and b4<hh and b5<hh:
   g=18
  elif b1>=hh and b2<hh and b3<hh and b4>=hh and b5<hh:
   g=19
  elif b1>=hh and b2<hh and b3<hh and b4<hh and b5>=hh:
   g=20
  elif b1<hh and b2>=hh and b3>=hh and b4<hh and b5<hh:
   g=21
  elif b1<hh and b2>=hh and b3<hh and b4>=hh and b5<hh:
   g=22
  elif b1<hh and b2>=hh and b3<hh and b4<hh and b5>=hh:
   g=23
  elif b1<hh and b2<hh and b3>=hh and b4>=hh and b5<hh:
   g=24
  elif b1<hh and b2<hh and b3>=hh and b4<hh and b5>=hh:
   g=25
  elif b1<hh and b2<hh and b3<hh and b4>=hh and b5>=hh:
   g=26
  elif b1>=hh and b2<hh and b3<hh and b4<hh and b5<hh:
   g=27
  elif b1<hh and b2>=hh and b3<hh and b4<hh and b5<hh:
   g=28
  elif b1<hh and b2<hh and b3>=hh and b4<hh and b5<hh:
   g=29
  elif b1<hh and b2<hh and b3<hh and b4>=hh and b5<hh:
   g=30
  elif b1<hh and b2<hh and b3<hh and b4<hh and b5>=hh:
   g=31
  else:
   c1=co2(a,49,cr,ob,bi)	#702.5
   c2=co2(a,51,cr,ob,bi)	#703.75
   c3=co2(a,53,cr,ob,bi)	#705
   c4=co2(a,55,cr,ob,bi)	#706.25
   c5=co2(a,57,cr,ob,bi)	#707.5
   if c1>=hh and c2>=hh and c3>=hh and c4>=hh and c5>=hh:
    g=32
   elif c1>=hh and c2>=hh and c3>=hh and c4>=hh and c5<hh:
    g=33
   elif c1>=hh and c2>=hh and c3>=hh and c4<hh and c5>=hh:
    g=34
   elif c1>=hh and c2>=hh and c3<hh and c4>=hh and c5>=hh:
    g=35
   elif c1>=hh and c2<hh and c3>=hh and c4>=hh and c5>=hh:
    g=36
   elif c1<hh and c2>=hh and c3>=hh and c4>=hh and c5>=hh:
    g=37
   elif c1>=hh and c2>=hh and c3>=hh and c4<hh and c5<hh:
    g=38
   elif c1>=hh and c2>=hh and c3<hh and c4>=hh and c5<hh:
    g=39
   elif c1>=hh and c2<hh and c3>=hh and c4>=hh and c5<hh:
    g=40
   elif c1<hh and c2>=hh and c3>=hh and c4>=hh and c5<hh:
    g=41
   elif c1>=hh and c2>=hh and c3<hh and c4<hh and c5>=hh:
    g=42
   elif c1>=hh and c2<hh and c3>=hh and c4<hh and c5>=hh:
    g=43
   elif c1<hh and c2>=hh and c3>=hh and c4<hh and c5>=hh:
    g=44
   elif c1>=hh and c2<hh and c3<hh and c4>=hh and c5>=hh:
    g=45
   elif c1<hh and c2>=hh and c3<hh and c4>=hh and c5>=hh:
    g=46
   elif c1<hh and c2<hh and c3>=hh and c4>=hh and c5>=hh:
    g=47
   elif c1>=hh and c2>=hh and c3<hh and c4<hh and c5<hh:
    g=48
   elif c1>=hh and c2<hh and c3>=hh and c4<hh and c5<hh:
    g=49
   elif c1>=hh and c2<hh and c3<hh and c4>=hh and c5<hh:
    g=50
   elif c1>=hh and c2<hh and c3<hh and c4<hh and c5>=hh:
    g=51
   elif c1<hh and c2>=hh and c3>=hh and c4<hh and c5<hh:
    g=52
   elif c1<hh and c2>=hh and c3<hh and c4>=hh and c5<hh:
    g=53
   elif c1<hh and c2>=hh and c3<hh and c4<hh and c5>=hh:
    g=54
   elif c1<hh and c2<hh and c3>=hh and c4>=hh and c5<hh:
    g=55
   elif c1<hh and c2<hh and c3>=hh and c4<hh and c5>=hh:
    g=56
   elif c1<hh and c2<hh and c3<hh and c4>=hh and c5>=hh:
    g=57
   elif c1>=hh and c2<hh and c3<hh and c4<hh and c5<hh:
    g=58
   elif c1<hh and c2>=hh and c3<hh and c4<hh and c5<hh:
    g=59
   elif c1<hh and c2<hh and c3>=hh and c4<hh and c5<hh:
    g=60
   elif c1<hh and c2<hh and c3<hh and c4>=hh and c5<hh:
    g=61
   elif c1<hh and c2<hh and c3<hh and c4<hh and c5>=hh:
    g=62
   else:
   d1=co2(a,65,cr,ob,bi)	#712.5
   d2=co2(a,67,cr,ob,bi)	#713.75
   d3=co2(a,69,cr,ob,bi)	#715
   d4=co2(a,71,cr,ob,bi)	#716.25
   d5=co2(a,73,cr,ob,bi)	#717.5
    if d1>=hh and d2>=hh and d3>=hh and d4>=hh and d5>=hh:
     g=63
    elif d1>=hh and d2>=hh and d3>=hh and d4>=hh and d5<hh:
     g=64
    elif d1>=hh and d2>=hh and d3>=hh and d4<hh and d5>=hh:
     g=65
    elif d1>=hh and d2>=hh and d3<hh and d4>=hh and d5>=hh:
     g=66
    elif d1>=hh and d2<hh and d3>=hh and d4>=hh and d5>=hh:
     g=67
    elif d1<hh and d2>=hh and d3>=hh and d4>=hh and d5>=hh:
     g=68
    elif d1>=hh and d2>=hh and d3>=hh and d4<hh and d5<hh:
     g=69
    elif d1>=hh and d2>=hh and d3<hh and d4>=hh and d5<hh:
     g=70
    elif d1>=hh and d2<hh and d3>=hh and d4>=hh and d5<hh:
     g=71
    elif d1<hh and d2>=hh and d3>=hh and d4>=hh and d5<hh:
     g=72
    elif d1>=hh and d2>=hh and d3<hh and d4<hh and d5>=hh:
     g=73
    elif d1>=hh and d2<hh and d3>=hh and d4<hh and d5>=hh:
     g=74
    elif d1<hh and d2>=hh and d3>=hh and d4<hh and d5>=hh:
     g=75
    elif d1>=hh and d2<hh and d3<hh and d4>=hh and d5>=hh:
     g=76
    elif d1<hh and d2>=hh and d3<hh and d4>=hh and d5>=hh:
     g=77
    elif d1<hh and d2<hh and d3>=hh and d4>=hh and d5>=hh:
     g=78
    elif d1>=hh and d2>=hh and d3<hh and d4<hh and d5<hh:
     g=79
    elif d1>=hh and d2<hh and d3>=hh and d4<hh and d5<hh:
     g=80
    elif d1>=hh and d2<hh and d3<hh and d4>=hh and d5<hh:
     g=81
    elif d1>=hh and d2<hh and d3<hh and d4<hh and d5>=hh:
     g=82
    elif d1<hh and d2>=hh and d3>=hh and d4<hh and d5<hh:
     g=83
    elif d1<hh and d2>=hh and d3<hh and d4>=hh and d5<hh:
     g=84
    elif d1<hh and d2>=hh and d3<hh and d4<hh and d5>=hh:
     g=85
    elif d1<hh and d2<hh and d3>=hh and d4>=hh and d5<hh:
     g=86
    elif d1<hh and d2<hh and d3>=hh and d4<hh and d5>=hh:
     g=87
    elif d1<hh and d2<hh and d3<hh and d4>=hh and d5>=hh:
     g=88
    elif d1>=hh and d2<hh and d3<hh and d4<hh and d5<hh:
     g=89
    elif d1<hh and d2>=hh and d3<hh and d4<hh and d5<hh:
     g=90
    elif d1<hh and d2<hh and d3>=hh and d4<hh and d5<hh:
     g=91
    elif d1<hh and d2<hh and d3<hh and d4>=hh and d5<hh:
     g=92
    elif d1<hh and d2<hh and d3<hh and d4<hh and d5>=hh:
     g=93
    else:
    e1=co2(a,94,cr,ob,bi)      #730.625
    e2=co2(a,96,cr,ob,bi)      #731.875
    e3=co2(a,98,cr,ob,bi)      #733.125
    e4=co2(a,100,cr,ob,bi)    #734.375
    e5=co2(a,102,cr,ob,bi)    #735.625
     if e1>=hh and e2>=hh and e3>=hh and e4>=hh and e5>=hh:
      g=94
     elif e1>=hh and e2>=hh and e3>=hh and e4>=hh and e5<hh:
      g=95
     elif e1>=hh and e2>=hh and e3>=hh and e4<hh and e5>=hh:
      g=96
     elif e1>=hh and e2>=hh and e3<hh and e4>=hh and e5>=hh:
      g=97
     elif e1>=hh and e2<hh and e3>=hh and e4>=hh and e5>=hh:
      g=98
     elif e1<hh and e2>=hh and e3>=hh and e4>=hh and e5>=hh:
      g=99
     elif e1>=hh and e2>=hh and e3>=hh and e4<hh and e5<hh:
      g=100
     elif e1>=hh and e2>=hh and e3<hh and e4>=hh and e5<hh:
      g=101
     elif e1>=hh and e2<hh and e3>=hh and e4>=hh and e5<hh:
      g=102
     elif e1<hh and e2>=hh and e3>=hh and e4>=hh and e5<hh:
      g=103
     elif e1>=hh and e2>=hh and e3<hh and e4<hh and e5>=hh:
      g=104
     elif e1>=hh and e2<hh and e3>=hh and e4<hh and e5>=hh:
      g=105
     elif e1<hh and e2>=hh and e3>=hh and e4<hh and e5>=hh:
      g=106
     elif e1>=hh and e2<hh and e3<hh and e4>=hh and e5>=hh:
      g=107
     elif e1<hh and e2>=hh and e3<hh and e4>=hh and e5>=hh:
      g=108
     elif e1<hh and e2<hh and e3>=hh and e4>=hh and e5>=hh:
      g=109
     elif e1>=hh and e2>=hh and e3<hh and e4<hh and e5<hh:
      g=110
     elif e1>=hh and e2<hh and e3>=hh and e4<hh and e5<hh:
      g=111
     elif e1>=hh and e2<hh and e3<hh and e4>=hh and e5<hh:
      g=112
     elif e1>=hh and e2<hh and e3<hh and e4<hh and e5>=hh:
      g=113
     elif e1<hh and e2>=hh and e3>=hh and e4<hh and e5<hh:
      g=114
     elif e1<hh and e2>=hh and e3<hh and e4>=hh and e5<hh:
      g=115
     elif e1<hh and e2>=hh and e3<hh and e4<hh and e5>=hh:
      g=116
     elif e1<hh and e2<hh and e3>=hh and e4>=hh and e5<hh:
      g=117
     elif e1<hh and e2<hh and e3>=hh and e4<hh and e5>=hh:
      g=118
     elif e1<hh and e2<hh and e3<hh and e4>=hh and e5>=hh:
      g=119
     elif e1>=hh and e2<hh and e3<hh and e4<hh and e5<hh:
      g=120
     elif e1<hh and e2>=hh and e3<hh and e4<hh and e5<hh:
      g=121
     elif e1<hh and e2<hh and e3>=hh and e4<hh and e5<hh:
      g=122
     elif e1<hh and e2<hh and e3<hh and e4>=hh and e5<hh:
      g=123
     elif e1<hh and e2<hh and e3<hh and e4<hh and e5>=hh:
      g=124
     else:
      g=125
  return(g)

 def pressure_co2(a,bi,cr,da,la,ob,pr,ta,te,to):
  b=level(a,cr,ob,bi,te)
  c=np.asarray(integrate_six(a,pr,te,da,la))
  if b==1:
   d1=find_nearest_index((integrate_five(a,28,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,30,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,32,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,34,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d5=find_nearest_index((integrate_five(a,36,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e4=c[d4]
   e5=c[d5]
   e6=[e1,e2,e3,e4,e5]
   e=trim_ave(e6)
  elif b==2:
   d1=find_nearest_index((integrate_five(a,28,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,30,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,32,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,34,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e4=c[d4]
   e6=[e1,e2,e3,e4]
   e=trim_ave(e6)
  elif b==3:
   d1=find_nearest_index((integrate_five(a,28,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,30,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,32,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,36,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e4=c[d4]
   e6=[e1,e2,e3,e4]
   e=trim_ave(e6)
  elif b==4:
   d1=find_nearest_index((integrate_five(a,28,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,30,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,34,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,36,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e4=c[d4]
   e6=[e1,e2,e3,e4]
   e=trim_ave(e6)
  elif b==5:
   d1=find_nearest_index((integrate_five(a,28,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,32,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,34,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,36,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e4=c[d4]
   e6=[e1,e2,e3,e4]
   e=trim_ave(e6)
  elif b==6:
   d1=find_nearest_index((integrate_five(a,30,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,32,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,34,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,36,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e4=c[d4]
   e6=[e1,e2,e3,e4]
   e=trim_ave(e6)
  elif b==7:
   d1=find_nearest_index((integrate_five(a,28,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,30,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,32,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==8:
   d1=find_nearest_index((integrate_five(a,28,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,30,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,34,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==9:
   d1=find_nearest_index((integrate_five(a,28,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,32,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,34,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==10:
   d1=find_nearest_index((integrate_five(a,30,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,32,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,34,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==11:
   d1=find_nearest_index((integrate_five(a,28,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,30,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,36,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==12:
   d1=find_nearest_index((integrate_five(a,28,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,32,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,36,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==13:
   d1=find_nearest_index((integrate_five(a,30,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,32,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,36,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==14:
   d1=find_nearest_index((integrate_five(a,28,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,34,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,36,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==15:
   d1=find_nearest_index((integrate_five(a,30,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,34,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,36,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==16:
   d1=find_nearest_index((integrate_five(a,32,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,34,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,36,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==17:
   d1=find_nearest_index((integrate_five(a,28,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,30,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e=(e1+e2)/2
  elif b==18:
   d1=find_nearest_index((integrate_five(a,28,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,32,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e3=c[d3]
   e=(e1+e3)/2
  elif b==19:
   d1=find_nearest_index((integrate_five(a,28,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,34,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e4=c[d4]
   e=(e1+e4)/2
  elif b==20:
   d1=find_nearest_index((integrate_five(a,28,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d5=find_nearest_index((integrate_five(a,36,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e5=c[d5]
   e=(e1+e5)/2
  elif b==21:
   d2=find_nearest_index((integrate_five(a,30,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,32,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e2=c[d2]
   e3=c[d3]
   e=(e2+e3)/2
  elif b==22:
   d2=find_nearest_index((integrate_five(a,30,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,34,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e2=c[d2]
   e4=c[d4]
   e=(e2+e4)/2
  elif b==23:
   d2=find_nearest_index((integrate_five(a,30,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d5=find_nearest_index((integrate_five(a,36,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e2=c[d2]
   e5=c[d5]
   e=(e2+e5)/2
  elif b==24:
   d3=find_nearest_index((integrate_five(a,32,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,34,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e3=c[d3]
   e4=c[d4]
   e=(e3+e4)/2
  elif b==25:
   d3=find_nearest_index((integrate_five(a,32,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d5=find_nearest_index((integrate_five(a,36,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e3=c[d3]
   e5=c[d5]
   e=(e3+e5)/2
  elif b==26:
   d4=find_nearest_index((integrate_five(a,34,53,ob,cr,ta,te,pr,da,la,bi)),0)
   d5=find_nearest_index((integrate_five(a,36,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e4=c[d4]
   e5=c[d5]
   e=(e4+e5)/2
  elif b==27:
   d1=find_nearest_index((integrate_five(a,28,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e=e1
  elif b==28:
   d2=find_nearest_index((integrate_five(a,30,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e2=c[d2]
   e=e2
  elif b==29:
   d3=find_nearest_index((integrate_five(a,32,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e3=c[d3]
   e=e3
  elif b==30:
   d4=find_nearest_index((integrate_five(a,34,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e4=c[d4]
   e=e4
  elif b==31:
   d5=find_nearest_index((integrate_five(a,36,53,ob,cr,ta,te,pr,da,la,bi)),0)
   e5=c[d5]
   e=e5
  elif b==32:
   d1=find_nearest_index((integrate_five(a,49,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,51,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,53,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,55,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d5=find_nearest_index((integrate_five(a,57,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e4=c[d4]
   e5=c[d5]
   e6=[e1,e2,e3,e4,e5]
   e=trim_ave(e6)
  elif b==33:
   d1=find_nearest_index((integrate_five(a,49,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,51,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,53,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,55,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e4=c[d4]
   e6=[e1,e2,e3,e4]
   e=trim_ave(e6)
  elif b==34:
   d1=find_nearest_index((integrate_five(a,49,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,51,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,53,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,57,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e4=c[d4]
   e6=[e1,e2,e3,e4]
   e=trim_ave(e6)
  elif b==35:
   d1=find_nearest_index((integrate_five(a,49,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,51,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,55,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,57,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e4=c[d4]
   e6=[e1,e2,e3,e4]
   e=trim_ave(e6)
  elif b==36:
   d1=find_nearest_index((integrate_five(a,49,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,53,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,55,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,57,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e4=c[d4]
   e6=[e1,e2,e3,e4]
   e=trim_ave(e6)
  elif b==37:
   d1=find_nearest_index((integrate_five(a,51,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,53,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,55,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,57,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e4=c[d4]
   e6=[e1,e2,e3,e4]
   e=trim_ave(e6)
  elif b==38:
   d1=find_nearest_index((integrate_five(a,49,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,51,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,53,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==39:
   d1=find_nearest_index((integrate_five(a,49,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,51,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,55,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==40:
   d1=find_nearest_index((integrate_five(a,49,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,53,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,55,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==41:
   d1=find_nearest_index((integrate_five(a,51,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,53,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,55,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==42:
   d1=find_nearest_index((integrate_five(a,49,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,51,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,57,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==43:
   d1=find_nearest_index((integrate_five(a,49,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,53,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,57,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==44:
   d1=find_nearest_index((integrate_five(a,51,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,53,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,57,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==45:
   d1=find_nearest_index((integrate_five(a,49,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,55,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,57,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==46:
   d1=find_nearest_index((integrate_five(a,51,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,55,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,57,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==47:
   d1=find_nearest_index((integrate_five(a,53,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,55,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,57,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==48:
   d1=find_nearest_index((integrate_five(a,49,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,51,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e=(e1+e2)/2
  elif b==49:
   d1=find_nearest_index((integrate_five(a,49,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,53,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e3=c[d3]
   e=(e1+e3)/2
  elif b==50:
   d1=find_nearest_index((integrate_five(a,49,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,55,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e4=c[d4]
   e=(e1+e4)/2
  elif b==51:
   d1=find_nearest_index((integrate_five(a,49,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d5=find_nearest_index((integrate_five(a,57,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e5=c[d5]
   e=(e1+e5)/2
  elif b==52:
   d2=find_nearest_index((integrate_five(a,51,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,53,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e2=c[d2]
   e3=c[d3]
   e=(e2+e3)/2
  elif b==53:
   d2=find_nearest_index((integrate_five(a,51,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,55,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e2=c[d2]
   e4=c[d4]
   e=(e2+e4)/2
  elif b==54:
   d2=find_nearest_index((integrate_five(a,51,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d5=find_nearest_index((integrate_five(a,57,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e2=c[d2]
   e5=c[d5]
   e=(e2+e5)/2
  elif b==55:
   d3=find_nearest_index((integrate_five(a,53,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,55,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e3=c[d3]
   e4=c[d4]
   e=(e3+e4)/2
  elif b==56:
   d3=find_nearest_index((integrate_five(a,53,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d5=find_nearest_index((integrate_five(a,57,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e3=c[d3]
   e5=c[d5]
   e=(e3+e5)/2
  elif b==57:
   d4=find_nearest_index((integrate_five(a,55,69,ob,cr,ta,te,pr,da,la,bi)),0)
   d5=find_nearest_index((integrate_five(a,57,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e4=c[d4]
   e5=c[d5]
   e=(e4+e5)/2
  elif b==58:
   d1=find_nearest_index((integrate_five(a,49,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e=e1
  elif b==59:
   d2=find_nearest_index((integrate_five(a,51,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e2=c[d2]
   e=e2
  elif b==60:
   d3=find_nearest_index((integrate_five(a,53,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e3=c[d3]
   e=e3
  elif b==61:
   d4=find_nearest_index((integrate_five(a,55,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e4=c[d4]
   e=e4
  elif b==62:
   d5=find_nearest_index((integrate_five(a,57,69,ob,cr,ta,te,pr,da,la,bi)),0)
   e5=c[d5]
   e=e5
  elif b==63:
   d1=find_nearest_index((integrate_five(a,65,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,67,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,69,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,71,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d5=find_nearest_index((integrate_five(a,73,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e4=c[d4]
   e5=c[d5]
   e6=[e1,e2,e3,e4,e5]
   e=trim_ave(e6)
  elif b==64:
   d1=find_nearest_index((integrate_five(a,65,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,67,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,69,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,71,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e4=c[d4]
   e6=[e1,e2,e3,e4]
   e=trim_ave(e6)
  elif b==65:
   d1=find_nearest_index((integrate_five(a,65,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,67,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,69,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,73,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e4=c[d4]
   e6=[e1,e2,e3,e4]
   e=trim_ave(e6)
  elif b==66:
   d1=find_nearest_index((integrate_five(a,65,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,67,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,71,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,73,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e4=c[d4]
   e6=[e1,e2,e3,e4]
   e=trim_ave(e6)
  elif b==67:
   d1=find_nearest_index((integrate_five(a,65,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,69,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,71,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,73,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e4=c[d4]
   e6=[e1,e2,e3,e4]
   e=trim_ave(e6)
  elif b==68:
   d1=find_nearest_index((integrate_five(a,67,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,69,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,71,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,73,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e4=c[d4]
   e6=[e1,e2,e3,e4]
   e=trim_ave(e6)
  elif b==69:
   d1=find_nearest_index((integrate_five(a,65,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,67,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,69,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==70:
   d1=find_nearest_index((integrate_five(a,65,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,67,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,71,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==71:
   d1=find_nearest_index((integrate_five(a,65,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,69,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,71,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==72:
   d1=find_nearest_index((integrate_five(a,67,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,69,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,71,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==73:
   d1=find_nearest_index((integrate_five(a,65,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,67,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,73,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==74:
   d1=find_nearest_index((integrate_five(a,65,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,69,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,73,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==75:
   d1=find_nearest_index((integrate_five(a,67,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,69,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,73,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==76:
   d1=find_nearest_index((integrate_five(a,65,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,71,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,73,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==77:
   d1=find_nearest_index((integrate_five(a,67,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,71,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,73,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==78:
   d1=find_nearest_index((integrate_five(a,69,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,71,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,73,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==79:
   d1=find_nearest_index((integrate_five(a,65,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,67,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e=(e1+e2)/2
  elif b==80:
   d1=find_nearest_index((integrate_five(a,65,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,69,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e3=c[d3]
   e=(e1+e3)/2
  elif b==81:
   d1=find_nearest_index((integrate_five(a,65,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,71,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e4=c[d4]
   e=(e1+e4)/2
  elif b==82:
   d1=find_nearest_index((integrate_five(a,65,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d5=find_nearest_index((integrate_five(a,73,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e5=c[d5]
   e=(e1+e5)/2
  elif b==83:
   d2=find_nearest_index((integrate_five(a,67,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,69,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e2=c[d2]
   e3=c[d3]
   e=(e2+e3)/2
  elif b==84:
   d2=find_nearest_index((integrate_five(a,67,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,71,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e2=c[d2]
   e4=c[d4]
   e=(e2+e4)/2
  elif b==85:
   d2=find_nearest_index((integrate_five(a,67,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d5=find_nearest_index((integrate_five(a,73,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e2=c[d2]
   e5=c[d5]
   e=(e2+e5)/2
  elif b==86:
   d3=find_nearest_index((integrate_five(a,69,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,71,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e3=c[d3]
   e4=c[d4]
   e=(e3+e4)/2
  elif b==87:
   d3=find_nearest_index((integrate_five(a,69,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d5=find_nearest_index((integrate_five(a,73,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e3=c[d3]
   e5=c[d5]
   e=(e3+e5)/2
  elif b==88:
   d4=find_nearest_index((integrate_five(a,71,98,ob,cr,ta,te,pr,da,la,bi)),0)
   d5=find_nearest_index((integrate_five(a,73,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e4=c[d4]
   e5=c[d5]
   e=(e4+e5)/2
  elif b==89:
   d1=find_nearest_index((integrate_five(a,65,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e=e1
  elif b==90:
   d2=find_nearest_index((integrate_five(a,67,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e2=c[d2]
   e=e2
  elif b==91:
   d3=find_nearest_index((integrate_five(a,69,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e3=c[d3]
   e=e3
  elif b==92:
   d4=find_nearest_index((integrate_five(a,71,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e4=c[d4]
   e=e4
  elif b==93:
   d5=find_nearest_index((integrate_five(a,73,98,ob,cr,ta,te,pr,da,la,bi)),0)
   e5=c[d5]
   e=e5
  elif b==94:
   d1=find_nearest_index((integrate_five(a,94,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,96,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,98,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,100,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d5=find_nearest_index((integrate_five(a,102,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e4=c[d4]
   e5=c[d5]
   e6=[e1,e2,e3+e4+e5]
   e=trim_ave(e6)
  elif b==95:
   d1=find_nearest_index((integrate_five(a,94,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,96,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,98,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,100,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e4=c[d4]
   e6=[e1,e2,e3,e4]
   e=trim_ave(e6)
  elif b==96:
   d1=find_nearest_index((integrate_five(a,94,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,96,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,98,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,102,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e4=c[d4]
   e6=[e1,e2,e3,e4]
   e=trim_ave(e6)
  elif b==97:
   d1=find_nearest_index((integrate_five(a,94,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,96,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,100,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,102,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e4=c[d4]
   e6=[e1,e2,e3,e4]
   e=trim_ave(e6)
  elif b==98:
   d1=find_nearest_index((integrate_five(a,94,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,98,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,100,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,102,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e4=c[d4]
   e6=[e1,e2,e3,e4]
   e=trim_ave(e6)
  elif b==99:
   d1=find_nearest_index((integrate_five(a,96,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,98,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,100,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,102,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e4=c[d4]
   e6=[e1,e2,e3,e4]
   e=trim_ave(e6)
  elif b==100:
   d1=find_nearest_index((integrate_five(a,94,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,96,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,98,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==101:
   d1=find_nearest_index((integrate_five(a,94,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,96,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,100,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==102:
   d1=find_nearest_index((integrate_five(a,94,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,98,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,100,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==103:
   d1=find_nearest_index((integrate_five(a,96,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,98,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,100,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d1]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==104:
   d1=find_nearest_index((integrate_five(a,94,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,96,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,102,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==105:
   d1=find_nearest_index((integrate_five(a,94,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,98,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,102,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==106:
   d1=find_nearest_index((integrate_five(a,96,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,98,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,102,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==107:
   d1=find_nearest_index((integrate_five(a,94,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,100,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,102,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==108:
   d1=find_nearest_index((integrate_five(a,96,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,100,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,102,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==109:
   d1=find_nearest_index((integrate_five(a,98,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,100,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,102,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e3=c[d3]
   e6=[e1,e2,e3]
   e=trim_ave(e6)
  elif b==110:
   d1=find_nearest_index((integrate_five(a,94,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d2=find_nearest_index((integrate_five(a,96,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e2=c[d2]
   e=(e1+e2)/2
  elif b==111:
   d1=find_nearest_index((integrate_five(a,94,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,98,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e3=c[d3]
   e=(e1+e3)/2
  elif b==112:
   d1=find_nearest_index((integrate_five(a,94,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,100,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e4=c[d4]
   e=(e1+e4)/2
  elif b==113:
   d1=find_nearest_index((integrate_five(a,94,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d5=find_nearest_index((integrate_five(a,102,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e5=c[d5]
   e=(e1+e5)/2
  elif b==114:
   d2=find_nearest_index((integrate_five(a,96,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d3=find_nearest_index((integrate_five(a,98,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e2=c[d2]
   e3=c[d3]
   e=(e2+e3)/2
  elif b==115:
   d2=find_nearest_index((integrate_five(a,96,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,100,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e2=c[d2]
   e4=c[d4]
   e=(e2+e4)/2
  elif b==116:
   d2=find_nearest_index((integrate_five(a,96,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d5=find_nearest_index((integrate_five(a,102,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e2=c[d2]
   e5=c[d5]
   e=(e2+e5)/2
  elif b==117:
   d3=find_nearest_index((integrate_five(a,98,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d4=find_nearest_index((integrate_five(a,100,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e3=c[d3]
   e4=c[d4]
   e=(e3+e4)/2
  elif b==118:
   d3=find_nearest_index((integrate_five(a,98,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d5=find_nearest_index((integrate_five(a,102,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e3=c[d3]
   e5=c[d5]
   e=(e3+e5)/2
  elif b==119:
   d4=find_nearest_index((integrate_five(a,100,122,ob,cr,ta,te,pr,da,la,bi)),0)
   d5=find_nearest_index((integrate_five(a,102,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e4=c[d4]
   e5=c[d5]
   e=(e4+e5)/2
  elif b==120:
   d1=find_nearest_index((integrate_five(a,94,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e1=c[d1]
   e=e1
  elif b==121:
   d2=find_nearest_index((integrate_five(a,96,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e2=c[d2]
   e=e2
  elif b==122:
   d3=find_nearest_index((integrate_five(a,98,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e3=c[d3]
   e=e3
  elif b==123:
   d4=find_nearest_index((integrate_five(a,100,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e4=c[d4]
   e=e4
  elif b==124:
   d5=find_nearest_index((integrate_five(a,102,122,ob,cr,ta,te,pr,da,la,bi)),0)
   e5=c[d5]
   e=e5
  else:
   e=-8888
  g1=ir_test_ctp(a,ob,te,pr)
  g2=integrate_one(a,216,ob,cr,bi)
  g3=to[a]
  g4=optics(a,ob,cr,te,bi,pr,g1)
  g5=optics(a,ob,cr,te,bi,pr,e)
  if g1>=50 and g1<=950 and e<=660 and e>=50 and g1>=e and g5<=1.3:
   h=e
  elif g1>=50 and g1<=950 and e<=660 and e>=50 and g1<e and g2>=0.5 and g4<=1.3:
   h=g1
  elif g1<50 and e>=50 and e<=660 and g5<=1.3:
   h=e
  elif g1>950 and e>=50 and e<=660 and g5<=1.3:
   h=e
  elif g1>=50 and g1<=950 and e<50 and g2>=0.5 and g4<=1.3:
   h=g1
  elif g1>=50 and g1<=950 and e>660 and g2>=0.5 and g4<=1.3:
   h=g1
  elif g1>950 and g2>=0.5:
   h=2100 #low cloud
  elif g2<=-0.33 and g3==1:
   h=2000 #clr certain
  elif g2<=0.5:
   h=2050 #clr probable
  else:
   h=-7777 #inconclusive
  if h==-7777:
   j=h
  elif h>1500:
   j=h
  else:
   h1=pr[a][:]
   h2=find_nearest_index(h1,h)+3
   if h2<=126:
    j=pr[a][h2]
   else:
    j=pr[a][126]
  return(j)

 #1=cloud but not used and not at nu_1 yet so "clear"
 #2=cloud but not used and after nu_1 so "cloudy"
 #3=cloud and used as nu_1
 #4=cloud and used as nu_2
 #5=cloud from ir test
 #6=clear
 #7=error code/inconclusive

 def chan_used(a,bi,cr,da,la,ob,pr,ta,te,to):
  b=level(a,cr,ob,bi,te)
  d=pressure_co2(a,bi,cr,da,la,ob,pr,ta,te,to)
  if b==1 and d>=0 and d<=660:
   c_28=3
   c_30=3
   c_32=3
   c_34=3
   c_36=3
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==2 and d>=0 and d<=660:
   c_28=3
   c_30=3
   c_32=3
   c_34=3
   c_36=1
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==3 and d>=0 and d<=660:
   c_28=3
   c_30=3
   c_32=3
   c_34=1
   c_36=3
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==4 and d>=0 and d<=660:
   c_28=3
   c_30=3
   c_32=1
   c_34=3
   c_36=3
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==5 and d>=0 and d<=660:
   c_28=3
   c_30=1
   c_32=3
   c_34=3
   c_36=3
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==6 and d>=0 and d<=660:
   c_28=1
   c_30=3
   c_32=3
   c_34=3
   c_36=3
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==7 and d>=0 and d<=660:
   c_28=3
   c_30=3
   c_32=3
   c_34=1
   c_36=1
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==8 and d>=0 and d<=660:
   c_28=3
   c_30=3
   c_32=1
   c_34=3
   c_36=1
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==9 and d>=0 and d<=660:
   c_28=3
   c_30=1
   c_32=3
   c_34=3
   c_36=1
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==10 and d>=0 and d<=660:
   c_28=1
   c_30=3
   c_32=3
   c_34=3
   c_36=1
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==11 and d>=0 and d<=660:
   c_28=3
   c_30=3
   c_32=1
   c_34=1
   c_36=3
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==12 and d>=0 and d<=660:
   c_28=3
   c_30=1
   c_32=3
   c_34=1
   c_36=3
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==13 and d>=0 and d<=660:
   c_28=1
   c_30=3
   c_32=3
   c_34=1
   c_36=3
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==14 and d>=0 and d<=660:
   c_28=3
   c_30=1
   c_32=1
   c_34=3
   c_36=3
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==15 and d>=0 and d<=660:
   c_28=1
   c_30=3
   c_32=1
   c_34=3
   c_36=3
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==16 and d>=0 and d<=660:
   c_28=1
   c_30=1
   c_32=3
   c_34=3
   c_36=3
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==17 and d>=0 and d<=660:
   c_28=3
   c_30=3
   c_32=1
   c_34=1
   c_36=1
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==18 and d>=0 and d<=660:
   c_28=3
   c_30=1
   c_32=3
   c_34=1
   c_36=1
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==19 and d>=0 and d<=660:
   c_28=3
   c_30=1
   c_32=1
   c_34=3
   c_36=1
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==20 and d>=0 and d<=660:
   c_28=3
   c_30=1
   c_32=1
   c_34=1
   c_36=3
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==21 and d>=0 and d<=660:
   c_28=1
   c_30=3
   c_32=3
   c_34=1
   c_36=1
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==22 and d>=0 and d<=660:
   c_28=1
   c_30=3
   c_32=1
   c_34=3
   c_36=1
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==23 and d>=0 and d<=660:
   c_28=1
   c_30=3
   c_32=1
   c_34=1
   c_36=3
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==24 and d>=0 and d<=660:
   c_28=1
   c_30=1
   c_32=3
   c_34=3
   c_36=1
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==25 and d>=0 and d<=660:
   c_28=1
   c_30=1
   c_32=3
   c_34=1
   c_36=3
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==26 and d>=0 and d<=660:
   c_28=1
   c_30=1
   c_32=1
   c_34=3
   c_36=3
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==27 and d>=0 and d<=660:
   c_28=3
   c_30=1
   c_32=1
   c_34=1
   c_36=1
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==28 and d>=0 and d<=660:
   c_28=1
   c_30=3
   c_32=1
   c_34=1
   c_36=1
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==29 and d>=0 and d<=660:
   c_28=1
   c_30=1
   c_32=3
   c_34=1
   c_36=1
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==30 and d>=0 and d<=660:
   c_28=1
   c_30=1
   c_32=1
   c_34=3
   c_36=1
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==31 and d>=0 and d<=660:
   c_28=1
   c_30=1
   c_32=1
   c_34=1
   c_36=3
   c_49,c_51,c_53,c_55,c_57=2,2,4,2,2
   c_65,c_67,c_69,c_71,c_73=2,2,2,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==32 and d>=0 and d<=660:
   c_49=3
   c_51=3
   c_53=3
   c_55=3
   c_57=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==33 and d>=0 and d<=660:
   c_49=3
   c_51=3
   c_53=3
   c_55=3
   c_57=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==34 and d>=0 and d<=660:
   c_49=3
   c_51=3
   c_53=3
   c_55=1
   c_57=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==35 and d>=0 and d<=660:
   c_49=3
   c_51=3
   c_53=1
   c_55=3
   c_57=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==36 and d>=0 and d<=660:
   c_49=3
   c_51=1
   c_53=3
   c_55=3
   c_57=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==37 and d>=0 and d<=660:
   c_49=1
   c_51=3
   c_53=3
   c_55=3
   c_57=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==38 and d>=0 and d<=660:
   c_49=3
   c_51=3
   c_53=3
   c_55=1
   c_57=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==39 and d>=0 and d<=660:
   c_49=3
   c_51=3
   c_53=1
   c_55=3
   c_57=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==40 and d>=0 and d<=660:
   c_49=3
   c_51=1
   c_53=3
   c_55=3
   c_57=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==41 and d>=0 and d<=660:
   c_49=1
   c_51=3
   c_53=3
   c_55=3
   c_57=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==42 and d>=0 and d<=660:
   c_49=3
   c_51=3
   c_53=1
   c_55=1
   c_57=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==43 and d>=0 and d<=660:
   c_49=3
   c_51=1
   c_53=3
   c_55=1
   c_57=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==44 and d>=0 and d<=660:
   c_49=1
   c_51=3
   c_53=3
   c_55=1
   c_57=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==45 and d>=0 and d<=660:
   c_49=3
   c_51=1
   c_53=1
   c_55=3
   c_57=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==46 and d>=0 and d<=660:
   c_49=1
   c_51=3
   c_53=1
   c_55=3
   c_57=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==47 and d>=0 and d<=660:
   c_49=1
   c_51=1
   c_53=3
   c_55=3
   c_57=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==48 and d>=0 and d<=660:
   c_49=3
   c_51=3
   c_53=1
   c_55=1
   c_57=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==49 and d>=0 and d<=660:
   c_49=3
   c_51=1
   c_53=3
   c_55=1
   c_57=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==50 and d>=0 and d<=660:
   c_49=3
   c_51=1
   c_53=1
   c_55=3
   c_57=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==51 and d>=0 and d<=660:
   c_49=3
   c_51=1
   c_53=1
   c_55=1
   c_57=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==52 and d>=0 and d<=660:
   c_49=1
   c_51=3
   c_53=3
   c_55=1
   c_57=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==53 and d>=0 and d<=660:
   c_49=1
   c_51=3
   c_53=1
   c_55=3
   c_57=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==54 and d>=0 and d<=660:
   c_49=1
   c_51=3
   c_53=1
   c_55=1
   c_57=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==55 and d>=0 and d<=660:
   c_49=1
   c_51=1
   c_53=3
   c_55=3
   c_57=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==56 and d>=0 and d<=660:
   c_49=1
   c_51=1
   c_53=3
   c_55=1
   c_57=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==57 and d>=0 and d<=660:
   c_49=1
   c_51=1
   c_53=1
   c_55=3
   c_57=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==58 and d>=0 and d<=660:
   c_49=3
   c_51=1
   c_53=1
   c_55=1
   c_57=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==59 and d>=0 and d<=660:
   c_49=1
   c_51=3
   c_53=1
   c_55=1
   c_57=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==60 and d>=0 and d<=660:
   c_49=1
   c_51=1
   c_53=3
   c_55=1
   c_57=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==61 and d>=0 and d<=660:
   c_49=1
   c_51=1
   c_53=1
   c_55=3
   c_57=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==62 and d>=0 and d<=660:
   c_49=1
   c_51=1
   c_53=1
   c_55=1
   c_57=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=2,2,4,2,2
   c_94,c_96,c_98,c_100,c_102=2,2,2,2,2
   c_122=2
  elif b==63 and d>=0 and d<=660:
   c_65=3
   c_67=3
   c_69=3
   c_71=3
   c_73=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==64 and d>=0 and d<=660:
   c_65=3
   c_67=3
   c_69=3
   c_71=3
   c_73=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==65 and d>=0 and d<=660:
   c_65=3
   c_67=3
   c_69=3
   c_71=1
   c_73=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==66 and d>=0 and d<=660:
   c_65=3
   c_67=3
   c_69=1
   c_71=3
   c_73=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==67 and d>=0 and d<=660:
   c_65=3
   c_67=1
   c_69=3
   c_71=3
   c_73=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==68 and d>=0 and d<=660:
   c_65=1
   c_67=3
   c_69=3
   c_71=3
   c_73=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==69 and d>=0 and d<=660:
   c_65=3
   c_67=3
   c_69=3
   c_71=1
   c_73=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==70 and d>=0 and d<=660:
   c_65=3
   c_67=3
   c_69=1
   c_71=3
   c_73=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==71 and d>=0 and d<=660:
   c_65=3
   c_67=1
   c_69=3
   c_71=3
   c_73=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==72 and d>=0 and d<=660:
   c_65=1
   c_67=3
   c_69=3
   c_71=3
   c_73=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==73 and d>=0 and d<=660:
   c_65=3
   c_67=3
   c_69=1
   c_71=1
   c_73=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==74 and d>=0 and d<=660:
   c_65=3
   c_67=1
   c_69=3
   c_71=1
   c_73=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==75 and d>=0 and d<=660:
   c_65=1
   c_67=3
   c_69=3
   c_71=1
   c_73=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==76 and d>=0 and d<=660:
   c_65=3
   c_67=1
   c_69=1
   c_71=3
   c_73=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==77 and d>=0 and d<=660:
   c_65=1
   c_67=3
   c_69=1
   c_71=3
   c_73=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==78 and d>=0 and d<=660:
   c_65=1
   c_67=1
   c_69=3
   c_71=3
   c_73=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==79 and d>=0 and d<=660:
   c_65=3
   c_67=3
   c_69=1
   c_71=1
   c_73=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==80 and d>=0 and d<=660:
   c_65=3
   c_67=1
   c_69=3
   c_71=1
   c_73=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==81 and d>=0 and d<=660:
   c_65=3
   c_67=1
   c_69=1
   c_71=3
   c_73=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==82 and d>=0 and d<=660:
   c_65=3
   c_67=1
   c_69=1
   c_71=1
   c_73=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==83 and d>=0 and d<=660:
   c_65=1
   c_67=3
   c_69=3
   c_71=1
   c_73=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==84 and d>=0 and d<=660:
   c_65=1
   c_67=3
   c_69=1
   c_71=3
   c_73=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==85 and d>=0 and d<=660:
   c_65=1
   c_67=3
   c_69=1
   c_71=1
   c_73=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==86 and d>=0 and d<=660:
   c_65=1
   c_67=1
   c_69=3
   c_71=3
   c_73=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==87 and d>=0 and d<=660:
   c_65=1
   c_67=1
   c_69=3
   c_71=1
   c_73=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==88 and d>=0 and d<=660:
   c_65=1
   c_67=1
   c_69=1
   c_71=3
   c_73=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==89 and d>=0 and d<=660:
   c_65=3
   c_67=1
   c_69=1
   c_71=1
   c_73=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==90 and d>=0 and d<=660:
   c_65=1
   c_67=3
   c_69=1
   c_71=1
   c_73=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==91 and d>=0 and d<=660:
   c_65=1
   c_67=1
   c_69=3
   c_71=1
   c_73=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==92 and d>=0 and d<=660:
   c_65=1
   c_67=1
   c_69=1
   c_71=3
   c_73=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==93 and d>=0 and d<=660:
   c_65=1
   c_67=1
   c_69=1
   c_71=1
   c_73=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_94,c_96,c_98,c_100,c_102=2,2,4,2,2
   c_122=2
  elif b==94 and d>=0 and d<=660:
   c_94=3
   c_96=3
   c_98=3
   c_100=3
   c_102=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==95 and d>=0 and d<=660:
   c_94=3
   c_96=3
   c_98=3
   c_100=3
   c_102=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==96 and d>=0 and d<=660:
   c_94=3
   c_96=3
   c_98=3
   c_100=1
   c_102=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==97 and d>=0 and d<=660:
   c_94=3
   c_96=3
   c_98=1
   c_100=3
   c_102=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==98 and d>=0 and d<=660:
   c_94=3
   c_96=1
   c_98=3
   c_100=3
   c_102=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==99 and d>=0 and d<=660:
   c_94=1
   c_96=3
   c_98=3
   c_100=3
   c_102=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==100 and d>=0 and d<=660:
   c_94=3
   c_96=3
   c_98=3
   c_100=1
   c_102=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==101 and d>=0 and d<=660:
   c_94=3
   c_96=3
   c_98=1
   c_100=3
   c_102=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==102 and d>=0 and d<=660:
   c_94=3
   c_96=1
   c_98=3
   c_100=3
   c_102=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==103 and d>=0 and d<=660:
   c_94=1
   c_96=3
   c_98=3
   c_100=3
   c_102=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==104 and d>=0 and d<=660:
   c_94=3
   c_96=3
   c_98=1
   c_100=1
   c_102=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==105 and d>=0 and d<=660:
   c_94=3
   c_96=1
   c_98=3
   c_100=1
   c_102=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==106 and d>=0 and d<=660:
   c_94=1
   c_96=3
   c_98=3
   c_100=1
   c_102=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==107 and d>=0 and d<=660:
   c_94=3
   c_96=1
   c_98=1
   c_100=3
   c_102=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==108 and d>=0 and d<=660:
   c_94=1
   c_96=3
   c_98=1
   c_100=3
   c_102=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==109 and d>=0 and d<=660:
   c_94=1
   c_96=1
   c_98=3
   c_100=3
   c_102=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==110 and d>=0 and d<=660:
   c_94=3
   c_96=3
   c_98=1
   c_100=1
   c_102=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==111 and d>=0 and d<=660:
   c_94=3
   c_96=1
   c_98=3
   c_100=1
   c_102=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==112 and d>=0 and d<=660:
   c_94=3
   c_96=1
   c_98=1
   c_100=3
   c_102=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==113 and d>=0 and d<=660:
   c_94=3
   c_96=1
   c_98=1
   c_100=1
   c_102=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==114 and d>=0 and d<=660:
   c_94=1
   c_96=3
   c_98=3
   c_100=1
   c_102=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==115 and d>=0 and d<=660:
   c_94=1
   c_96=3
   c_98=1
   c_100=3
   c_102=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==116 and d>=0 and d<=660:
   c_94=1
   c_96=3
   c_98=1
   c_100=1
   c_102=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==117 and d>=0 and d<=660:
   c_94=1
   c_96=1
   c_98=3
   c_100=3
   c_102=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==118 and d>=0 and d<=660:
   c_94=1
   c_96=1
   c_98=3
   c_100=1
   c_102=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==119 and d>=0 and d<=660:
   c_94=1
   c_96=1
   c_98=1
   c_100=3
   c_102=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==120 and d>=0 and d<=660:
   c_94=3
   c_96=1
   c_98=1
   c_100=1
   c_102=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==121 and d>=0 and d<=660:
   c_94=1
   c_96=3
   c_98=1
   c_100=1
   c_102=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==122 and d>=0 and d<=660:
   c_94=1
   c_96=1
   c_98=3
   c_100=1
   c_102=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==123 and d>=0 and d<=660:
   c_94=1
   c_96=1
   c_98=1
   c_100=3
   c_102=1
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif b==124 and d>=0 and d<=660:
   c_94=1
   c_96=1
   c_98=1
   c_100=1
   c_102=3
   c_28,c_30,c_32,c_34,c_36=1,1,1,1,1
   c_49,c_51,c_53,c_55,c_57=1,1,1,1,1
   c_65,c_67,c_69,c_71,c_73=1,1,1,1,1
   c_122=4
  elif d>=0 and d<=1000:
   c_28,c_30,c_32,c_34,c_36=5,5,5,5,5
   c_49,c_51,c_53,c_55,c_57=5,5,5,5,5
   c_65,c_67,c_69,c_71,c_73=5,5,5,5,5
   c_94,c_96,c_98,c_100,c_102=5,5,5,5,5
   c_122=5
  elif d>=2000:
   c_28,c_30,c_32,c_34,c_36=6,6,6,6,6
   c_49,c_51,c_53,c_55,c_57=6,6,6,6,6
   c_65,c_67,c_69,c_71,c_73=6,6,6,6,6
   c_94,c_96,c_98,c_100,c_102=6,6,6,6,6
   c_122=6
  else:
   c_28,c_30,c_32,c_34,c_36=7,7,7,7,7
   c_49,c_51,c_53,c_55,c_57=7,7,7,7,7
   c_65,c_67,c_69,c_71,c_73=7,7,7,7,7
   c_94,c_96,c_98,c_100,c_102=7,7,7,7,7
   c_122=7
  return([c_28,c_30,c_32,c_34,c_36,c_49,c_51,c_53,c_55,c_57,c_65,c_67,c_69,c_71,c_73,c_94,c_96,c_98,c_100,c_102,c_122])

 b1=int(st1)
 b2=int(st2)
 b3=int(sp1)
 b4=int(sp2)
 b5=int(dt1)
 b6=int(dt2)
 b7=int(mon1)
 b8=int(ye1)
 if b1>=0 and b1<=9:
  c1="0"+str(b1)
 else:
  c1=str(b1)
 if b2==0:
  c2="00"
 else:
  c2="30"
 if b3>=0 and b3<=9:
  c3="0"+str(b3)
 else:
  c3=str(b3)
 if b4==0:
  c4="00"
 else:
  c4="30"
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
 st=c1+c2
 sp=c3+c4
 dt=c5
 dp=c6
 mon=c7
 yer=c8
 f1=cdf.Dataset("/apollo/jung/snebuda/brianne/2021/cris.431.n20."+yer+mon+dt+"_"+st+"00."+yer+mon+dp+"_"+sp+"00.bt.nc")
 f2=cdf.Dataset("/apollo/jung/snebuda/brianne/2021/cris.431.n20."+yer+mon+dt+"_"+st+"00."+yer+mon+dp+"_"+sp+"00.wf.nc")
 f3=cdf.Dataset("/apollo/jung/snebuda/brianne/2021/cris.431.n20."+yer+mon+dt+"_"+st+"00."+yer+mon+dp+"_"+sp+"00.state.nc")
 f4=cdf.Dataset("/apollo/jung/snebuda/brianne/2021/cris.431.n20."+yer+mon+dt+"_"+st+"00."+yer+mon+dp+"_"+sp+"00.te.nc")
 day_time=(yer+"_"+mon+"_"+dt+"_"+st+"_"+sp)
 output=("/apollo/jung/bandersen/co2_slicing_full/"+yer+"_"+mon+"_"+dt+"/co2_slicing_output_"+day_time+".txt")
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
 si=np.size(la)
 viirs_ctp=[]
 viirs_cot=[]
 co2_ctp=[]
 co2_cot=[]
 co2_level=[]
 trop_all=[]
 ir_test=[]
 use_28=[]
 use_30=[]
 use_32=[]
 use_34=[]
 use_36=[]
 use_49=[]
 use_51=[]
 use_53=[]
 use_55=[]
 use_57=[]
 use_65=[]
 use_67=[]
 use_69=[]
 use_71=[]
 use_73=[]
 use_94=[]
 use_96=[]
 use_98=[]
 use_100=[]
 use_102=[]
 use_122=[]
 obs_crtm_28=[]
 obs_crtm_30=[]
 obs_crtm_32=[]
 obs_crtm_34=[]
 obs_crtm_36=[]
 obs_crtm_49=[]
 obs_crtm_51=[]
 obs_crtm_53=[]
 obs_crtm_55=[]
 obs_crtm_57=[]
 obs_crtm_65=[]
 obs_crtm_67=[]
 obs_crtm_69=[]
 obs_crtm_71=[]
 obs_crtm_73=[]
 obs_crtm_94=[]
 obs_crtm_96=[]
 obs_crtm_98=[]
 obs_crtm_100=[]
 obs_crtm_102=[]
 obs_crtm_122=[]
 obs_crtm_216=[]
 obs_28=[]
 obs_30=[]
 obs_32=[]
 obs_34=[]
 obs_36=[]
 obs_49=[]
 obs_51=[]
 obs_53=[]
 obs_55=[]
 obs_57=[]
 obs_65=[]
 obs_67=[]
 obs_69=[]
 obs_71=[]
 obs_73=[]
 obs_94=[]
 obs_96=[]
 obs_98=[]
 obs_100=[]
 obs_102=[]
 obs_122=[]
 obs_216=[]
 for a in range(0,si):
  b=h2p(pv[a])
  c=tc[a]
  d=pressure_co2(a,bi,cr,da,la,ob,pr,ta,te,to)
# e=optics(a,ob,cr,te,bi)
  e=optics(a,ob,cr,te,bi,pr,d)
  g=level(a,cr,ob,bi,te)
  h=pr[a][trop(a,pr,te,da,la)]
  j=ir_test_ctp(a,ob,te,pr)
  k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20,k21=chan_used(a,bi,cr,da,la,ob,pr,ta,te,to)
  l1=co2(a,28,cr,ob,bi)
  l2=co2(a,30,cr,ob,bi)
  l3=co2(a,32,cr,ob,bi)
  l4=co2(a,34,cr,ob,bi)
  l5=co2(a,36,cr,ob,bi)
  l6=co2(a,49,cr,ob,bi)
  l7=co2(a,51,cr,ob,bi)
  l8=co2(a,53,cr,ob,bi)
  l9=co2(a,55,cr,ob,bi)
  l10=co2(a,57,cr,ob,bi)
  l11=co2(a,65,cr,ob,bi)
  l12=co2(a,67,cr,ob,bi)
  l13=co2(a,69,cr,ob,bi)
  l14=co2(a,71,cr,ob,bi)
  l15=co2(a,73,cr,ob,bi)
  l16=co2(a,94,cr,ob,bi)
  l17=co2(a,96,cr,ob,bi)
  l18=co2(a,98,cr,ob,bi)
  l19=co2(a,100,cr,ob,bi)
  l20=co2(a,102,cr,ob,bi)
  l21=co2(a,122,cr,ob,bi)
  l22=co2(a,216,cr,ob,bi)
  m1=ob[a][28]
  m2=ob[a][30]
  m3=ob[a][32]
  m4=ob[a][34]
  m5=ob[a][36]
  m6=ob[a][49]
  m7=ob[a][51]
  m8=ob[a][53]
  m9=ob[a][55]
  m10=ob[a][57]
  m11=ob[a][65]
  m12=ob[a][67]
  m13=ob[a][69]
  m14=ob[a][71]
  m15=ob[a][73]
  m16=ob[a][94]
  m17=ob[a][96]
  m18=ob[a][98]
  m19=ob[a][100]
  m20=ob[a][102]
  m21=ob[a][122]
  m22=ob[a][216]
  viirs_ctp.append(b)
  viirs_cot.append(c)
  co2_ctp.append(d)
  co2_cot.append(e)
  co2_level.append(g)
  trop_all.append(h)
  ir_test.append(j)
  use_28.append(k1)
  use_30.append(k2)
  use_32.append(k3)
  use_34.append(k4)
  use_36.append(k5)
  use_49.append(k6)
  use_51.append(k7)
  use_53.append(k8)
  use_55.append(k9)
  use_57.append(k10)
  use_65.append(k11)
  use_67.append(k12)
  use_69.append(k13)
  use_71.append(k14)
  use_73.append(k15)
  use_94.append(k16)
  use_96.append(k17)
  use_98.append(k18)
  use_100.append(k19)
  use_102.append(k20)
  use_122.append(k21)
  obs_crtm_28.append(l1)
  obs_crtm_30.append(l2)
  obs_crtm_32.append(l3)
  obs_crtm_34.append(l4)
  obs_crtm_36.append(l5)
  obs_crtm_49.append(l6)
  obs_crtm_51.append(l7)
  obs_crtm_53.append(l8)
  obs_crtm_55.append(l9)
  obs_crtm_57.append(l10)
  obs_crtm_65.append(l11)
  obs_crtm_67.append(l12)
  obs_crtm_69.append(l13)
  obs_crtm_71.append(l14)
  obs_crtm_73.append(l15)
  obs_crtm_94.append(l16)
  obs_crtm_96.append(l17)
  obs_crtm_98.append(l18)
  obs_crtm_100.append(l19)
  obs_crtm_102.append(l20)
  obs_crtm_122.append(l21)
  obs_crtm_216.append(l22)
  obs_28.append(m1)
  obs_30.append(m2)
  obs_32.append(m3)
  obs_34.append(m4)
  obs_36.append(m5)
  obs_49.append(m6)
  obs_51.append(m7)
  obs_53.append(m8)
  obs_55.append(m9)
  obs_57.append(m10)
  obs_65.append(m11)
  obs_67.append(m12)
  obs_69.append(m13)
  obs_71.append(m14)
  obs_73.append(m15)
  obs_94.append(m16)
  obs_96.append(m17)
  obs_98.append(m18)
  obs_100.append(m19)
  obs_102.append(m20)
  obs_122.append(m21)
  obs_216.append(m22)
 file=open(output,"w+")
  with open(output,'wb') as f:
   for a in range(0,si):
    z=np.array([la[a],lo[a],to[a],fo[a],fv[a],viirs_ctp[a],viirs_cot[a],co2_ctp[a],co2_cot[a],co2_level[a],trop_all[a],ir_test[a],use_28[a],use_30[a],use_32[a],use_34[a],use_36[a],use_49[a],use_51[a],use_53[a],use_55[a],use_57[a],use_65[a],use_67[a],use_69[a],use_71[a],use_73[a],use_94[a],use_96[a],use_98[a],use_100[a],use_102[a],use_122[a],obs_crtm_28[a],obs_crtm_30[a],obs_crtm_32[a],obs_crtm_34[a],obs_crtm_36[a],obs_crtm_49[a],obs_crtm_51[a],obs_crtm_53[a],obs_crtm_55[a],obs_crtm_57[a],obs_crtm_65[a],obs_crtm_67[a],obs_crtm_69[a],obs_crtm_71[a],obs_crtm_73[a],obs_crtm_94[a],obs_crtm_96[a],obs_crtm_98[a],obs_crtm_100[a],obs_crtm_102[a],obs_crtm_122[a],obs_crtm_216[a],obs_28[a],obs_30[a],obs_32[a],obs_34[a],obs_36[a],obs_49[a],obs_51[a],obs_53[a],obs_55[a],obs_57[a],obs_65[a],obs_67[a],obs_69[a],obs_71[a],obs_73[a],obs_94[a],obs_96[a],obs_98[a],obs_100[a],obs_102[a],obs_122[a],obs_216[a]])
    np.savetxt(f,z,fmt=' 1.7f',newline=" ")
    np.savetxt(f,[-8888],fmt=' d',newline=" ")
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
 return("doneer")
