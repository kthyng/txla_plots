'''
Plot mean sea surface height in the model domain from the model output. 
Made for the shelf transport paper.
'''

import numpy as np
import matplotlib.pyplot as plt
import glob
import netCDF4 as netCDF
import tracpy
import tracpy.plotting
import matplotlib as mpl
import pdb
import op
from matplotlib import ticker
from matplotlib.mlab import find
import os


mpl.rcParams.update({'font.size': 14})
mpl.rcParams['font.sans-serif'] = 'Arev Sans, Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Helvetica, Avant Garde, sans-serif'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.cal'] = 'cursive'
mpl.rcParams['mathtext.rm'] = 'sans'
mpl.rcParams['mathtext.tt'] = 'monospace'
mpl.rcParams['mathtext.it'] = 'sans:italic'
mpl.rcParams['mathtext.bf'] = 'sans:bold'
mpl.rcParams['mathtext.sf'] = 'sans'
mpl.rcParams['mathtext.fallback_to_cm'] = 'True'


whichtime = 'interannual' # 'seasonal' or 'interannual'
whicharea = 'summer' # 'winter' or 'summer'

loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
grid = tracpy.inout.readgrid(loc, usebasemap=True)

d = netCDF.Dataset(loc)
t = d.variables['ocean_time'][:]
units = d.variables['ocean_time'].units
dates = netCDF.num2date(t, units)
years = np.asarray([dates[i].year for i in xrange(len(dates))])
months = np.asarray([dates[i].month for i in xrange(len(dates))])

# desired time indices for average, 2004-2011 and Jan or Feb
tinds = (years>=2004) * (years<=2011) * ((months==1) + (months==2))
zetaw = d.variables['zeta'][tinds,:,:]
# zetaw.set_fill_value(np.nan)
# zetaw = zetaw.filled()
# zmaxw = np.nanmax(abs(zetaw))
tinds = (years>=2004) * (years<=2011) * ((months==7) + (months==8))
zetas = d.variables['zeta'][tinds,:,:]
# zetas.set_fill_value(np.nan)
# zetas = zetas.filled()
# zmaxs = np.nanmax(abs(zetas))

# levels = np.arange(-0.9, 1.1, 0.2)
# levels = np.arange(-1.8, 2.2, 0.4)
# # mean
# levels = np.arange(-.18, .22, 0.04)
# abs mean
# levels = np.arange(0, .22, 0.02)
# # min mean
# levels = np.arange(-.22, 0, 0.02)
# mean
levels = np.arange(0, .22, 0.02)

sshfname = 'calcs/ssh.npz'

fig, axarr = plt.subplots(1,2)
fig.set_size_inches(13.675, 6.6125)
fig.subplots_adjust(left=0.04, bottom=0.15, right=1.0, top=0.96, wspace=0.07, hspace=0.04)

for i, ax in enumerate(axarr):

    # Titles for subplots
    if i==0:

        tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2))
        ax.set_title('Winter')
        # ind = zetaw<0
        # ax.contourf(grid['xr'].T, grid['yr'].T, zetaw[ind].mean(axis=0), cmap='Blues', levels=levels, extend='min')
        # ax.contourf(grid['xr'].T, grid['yr'].T, abs(zetaw).mean(axis=0), cmap='Reds', levels=levels, extend='max')
        # ax.contourf(grid['xr'].T, grid['yr'].T, zetaw.mean(axis=0), cmap='RdBu_r', levels=levels, extend='both')
        ax.contourf(grid['xr'].T, grid['yr'].T, zetaw.mean(axis=0)+0.18, cmap='Reds', levels=levels, extend='both')

    elif i==1:

        tracpy.plotting.background(grid=grid, ax=ax, parslabels=[0,0,0,0], mers=np.arange(-100, -80, 2))
        ax.set_title('Summer')
        # ind = zetas<0
        # mappable = ax.contourf(grid['xr'].T, grid['yr'].T, zetas[ind].mean(axis=0), cmap='Reds', levels=levels, extend='min')
        # mappable = ax.contourf(grid['xr'].T, grid['yr'].T, abs(zetas).mean(axis=0), cmap='Reds', levels=levels, extend='max')
        # mappable = ax.contourf(grid['xr'].T, grid['yr'].T, zetas.mean(axis=0), cmap='RdBu_r', levels=levels, extend='both')
        mappable = ax.contourf(grid['xr'].T, grid['yr'].T, zetas.mean(axis=0), cmap='Reds', levels=levels, extend='both')

    # zmax = np.nanmax((zmaxw,zmaxs))

# np.savez(sshfname, zetaw=zetaw.filled(), zetas=zetas.filled(), xr=grid['xr'].T, yr=grid['yr'].T)

# Horizontal colorbar below plot
cax = fig.add_axes([0.25, 0.075, 0.5, 0.02]) #colorbar axes
cb = plt.colorbar(mappable, cax=cax, orientation='horizontal')
cb.set_label('Mean sea surface height [m]')

# fig.savefig('figures/ssh/seasonal_min-mean.png', bbox_inches='tight')
# fig.savefig('figures/ssh/seasonal_abs-mean.png', bbox_inches='tight')
# fig.savefig('figures/ssh/seasonal_mean.png', bbox_inches='tight')
fig.savefig('figures/ssh/seasonal_mean-shifted.png', bbox_inches='tight')


