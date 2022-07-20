#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 15:38:24 2022

@author: francescalappin


Takes a netcdf of saved profile objects 
"""

from datetime import datetime, timedelta

import cmocean
import matplotlib.pyplot as plt
import scipy.interpolate as sc
import numpy as np
from metpy import calc
from metpy.units import units
import matplotlib.dates as mdates
from netCDF4 import Dataset
import cartopy.crs as ccrs
import netCDF4 as nc
from glob import glob
import metpy
import matplotlib.dates as datenum
import xarray as xr
from matplotlib import cm, colors
import pandas as pd

vars = {'theta': ["Potential Temperature", 'theta', 'K', cmocean.cm.thermal,
                  1.0],
        'temp': ["Temperature", 'temp', 'K', cmocean.cm.thermal, 1.0],
        'T_d': ["Dewpoint Temperature", 'T_d', '$^\circ$C', cmocean.cm.haline,
                1.0],
        'dewp': ["Dewpoint Temperature", 'T_d', '$^\circ$C', cmocean.cm.haline,
                 1.0],
        'r': ["Mixing Ratio", 'mixing_ratio', 'g Kg$^{-1}$', cmocean.cm.haline,
              0.5],
        'mr': ["Mixing Ratio", 'mixing_ratio', 'Kg Kg$^{-1}$', cmocean.cm.haline_r,
               0.5],
        'q': ["Specific Humidity", 'q', 'g Kg$^{-1}$', cmocean.cm.haline, 0.5],
        'rh': ["Relative Humidity", 'rh', '%', cmocean.cm.haline, 5.0],
        'speed': ["Wind Speed", 'speed', 'm s$^{-1}$', cmocean.cm.speed, 5.0],
        'ws': ["Wind Speed", 'speed', 'm s$^{-1}$', 'gist_stern_r', 5.0],
        'u': ["U", 'u', 'm s$^{-1}$', cmocean.cm.speed, 5.0],
        'v': ["V", 'v', 'm s$^{-1}$', cmocean.cm.speed, 5.0],
        'dir': ["Wind Direction", 'dir', '$^\circ$', cmocean.cm.phase, 360.,
                'wind'],
        'pres': ["Pressure", 'pres', 'Pa', cmocean.cm.haline, 15.0],
        'p': ["Pressure", 'pres', 'Pa', cmocean.cm.haline, 15.0],
        'alt': ["Altitude", 'alt', 'm', cmocean.cm.haline, 10.0]}

file = '/Users/francesca.lappin/CASS/TRACER-AQ/UAS-Profiles/processed/McCorN937UA.c1.20220621.cdf' 
#%%  
#Read in data and set to standard grid
var_l = ['temp','mr','speed','dir']   #order in how you want plotted
  
data_dic = {}  # also unitless
times_t=[]
z = []
for var_i in var_l:
    data_dic[var_i] = []

all_profs = Dataset(file,'r')
prof_names = list(all_profs.groups.keys())

for k in prof_names:
    # Get data from Profile objects
    nc_dat = all_profs[k]
    if (nc_dat['alt'][-1] < 300):
        continue                #This should skip short profiles (good for tracer but not for perils)
        
    t0 = nc_dat['time'][:]    
    t = list(nc.num2date(t0, units= 'microseconds since 2010-01-01 00:00:00:00'))
    times_t.append(t)
    
    z.append(nc_dat["alt"][:])
    for var_i in var_l:
        data_dic[var_i].append(list(nc_dat[var_i][:]))

max_len = 0
which_i = -1
for i in range(len(z)):             
    if len(z[i]) > max_len:
        max_len = len(z[i])
        which_i = i
          
z_t = z[which_i]

time_flat = np.array(times_t[0])
for p in range(len(times_t)):
    if p > 0:
        time_flat = np.concatenate((time_flat, times_t[p]))
timerange = datenum.drange(np.nanmin(time_flat),
                           np.nanmax(time_flat),
                           (np.nanmax(time_flat)
                            -np.nanmin(time_flat))/100)
for i in range(len(times_t)):
    diff = max_len - len(times_t[i])
    for j in range(diff):
        times_t[i].append(None)
        for var_i in var_l:
            data_dic[var_i][i].append(np.nan)
time = times_t


for p in range(len(time)):
    for i in range(len(time[p])):
        try:
            time[p][i] = datenum.date2num(time[p][i])
        except AttributeError:
            time[p][i] = np.nan

# Switch to arrays to make indexing easier
time = np.array(time, dtype=float)
z = np.array(z_t)
XX, YY = np.meshgrid(timerange, z)

# Prepare for interpolation
data_grid = {}
fig = None
for var_i in var_l:
    data_grid[var_i] = np.full_like(XX, np.nan)
    data_dic[var_i] = np.array(data_dic[var_i])


for var_i in var_l:
        for i in range(len(z)):
            a = list(np.array(time[:,i]).ravel())
            j = 0
            while j < len(a):
                if np.isnan((a[j])):
                    a.remove(a[j])
                else:
                    j += 1
            if len(a) < 2:
                continue
            interp_fun = sc.interp1d(np.array(time[:,i]).ravel(), data_dic[var_i][:,i],
                                  fill_value='extrapolate', kind='linear')      #cubic interp was going too far beyond data
            data_grid[var_i][i, :] = interp_fun(XX[i, :])    
        
#%%    
num_figs = len(var_l)       #plz use even numbers
num_rows = int(num_figs/2)

fig, ax = plt.subplots(nrows=num_rows, ncols=2, figsize=(15,17))

prof_start = np.array([e[0] for e in times_t])
xx, yy = np.meshgrid(prof_start, z)         #to compare against data_dic un-interp data    

plt.ylim((z[0], z[-1]))       

for i, ax in enumerate(fig.axes):
    cfax = ax.pcolormesh(XX, YY,data_grid[var_l[i]],
                         cmap=vars[var_l[i]][3])
    cbar = plt.colorbar(cfax, ax=ax, pad=0.01, )     
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45)
    ax.xaxis.set_major_formatter(datenum.DateFormatter('%H:%M:%S'))
    cbar.set_label(vars[var_l[i]][0] + " (" + str(vars[var_l[i]][2]) + ")",
                   rotation=270, fontsize=20, labelpad=30)
    for p_times in times_t:
        ax.scatter(p_times, z, c='black', s=0.5)
 
plt.suptitle(time_flat[0].strftime('%Y/%m/%d'),x=0.04,y=0.985,fontsize=24,ha='left')
plt.tight_layout()
#fig.savefig(('/Users/francesca.lappin/CASS/'+ time_flat[0].strftime('%Y%m%d') + 'cop4panel.png'))
    
#%%

pressure = data_dic['pres'] / 100           #lets get it in hpa
density = calc.density(pressure* units.mbar,data_dic['temp']* units.kelvin,data_dic['mr'])
fig,ax = plt.subplots(1,1,figsize=(12,8))
prs_anom_lev = np.arange(-0.5,0.5,0.05)
prs_lev = np.arange(950,1017,3)
dens_lev = np.arange(1.05,1.18,0.005)
cd=ax.contourf(xx.T,yy.T,np.gradient(pressure,axis=0),cmap='bwr',levels=prs_anom_lev,extend='both')
cc=ax.contour(xx.T,yy.T,density.magnitude, colors='k',levels=dens_lev,extend='both')
#cc=ax.contourf(XX,YY,density.magnitude, cmap=‘magma_r’,levels=dens_lev,extend=‘both’)
ax.clabel(cc,dens_lev, inline=True, fontsize=10)
cbar=plt.colorbar(cd)
cbar.set_label('Pressure Gradient [hPa]', size=16)
ax.set_ylim(0,600)
ax.set_xlabel('Time (UTC)')
ax.set_ylabel('Height AGL (km)')
ax.xaxis.set_major_locator(mdates.HourLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
plt.setp(ax.xaxis.get_majorticklabels(), rotation=45)
plt.title(f"Density and p’ {time_flat[0]:%Y-%m-%d}", loc='left',fontsize=16)




            
            
            
            