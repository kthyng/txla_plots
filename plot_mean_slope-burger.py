'''
Plot mean slope burger number in the model domain from the model output. 
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
import gsw
import octant
from scipy import ndimage


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

g = 9.81
rho0 = 1023. # kg/m^3

whichtime = 'seasonal' # 'seasonal' or 'interannual'
whicharea = 'summer' # 'winter' or 'summer'

loc = 'http://barataria.tamu.edu:6060/thredds/dodsC/NcML/txla_nesting6.nc'
grid = tracpy.inout.readgrid(loc, usebasemap=True)

d = netCDF.Dataset(loc)
t = d.variables['ocean_time'][:]
units = d.variables['ocean_time'].units
dates = netCDF.num2date(t, units)
years = np.asarray([dates[i].year for i in xrange(len(dates))])
months = np.asarray([dates[i].month for i in xrange(len(dates))])

fname1 = 'calcs/f-alpha.npz'
fname2 = 'calcs/' + whichtime + 'N.npz'

# Things that don't change in time
f = d.variables['f'][:] # Coriolis
pm = d.variables['pm'][:]; pn = d.variables['pn'][:]
h = d.variables['h'][:]
h_smoothed = ndimage.filters.gaussian_filter(h, 20)
dhdy, dhdx = np.gradient(h_smoothed)
alpha = np.sqrt((dhdx*pm)**2 + (dhdy*pn)**2) # dimensional bottom slope, across-shelf
np.savez(fname1, f=f, alpha=alpha, xr=grid['xr'].T, yr=grid['yr'].T)

Years = np.arange(2004,2012)

if whichtime == 'interannual':

    N = np.zeros((Years.size,2,f.shape)) # year x season x grid
    if not os.path.exists(fname2):
        for i,Year in enumerate(Years): # loop through years
            # Winter
            tinds =  find((Year==years) * ((months==1) + (months==2))) # loop through season
            for tind in tinds:
                salt = d.variables['salt'][tind,:,:,:]
                temp = d.variables['temp'][tind,:,:,:]
                rho = gsw.rho(salt, temp, 0)
                zeta = d.variables['zeta'][tind,:,:]
                zwt = octant.depths.get_zw(d.variables['Vtransform'][:][0], d.variables['Vstretching'][:][0], 
                            salt.shape[0], d.variables['theta_s'][:][0], d.variables['theta_b'][:][0], 
                                h, d.variables['hc'][:][0], zeta=zeta, Hscale=3)
                N[i,0,:,:] = N[i,0,:,:] + np.ma.median(np.sqrt(-g/rho0 * ((rho[2:,:,:]-rho[:-2,:,:])/(zwt[2:,:,:]-zwt[:-2,:,:]))), axis=0) # just save median
            N[i,0,:,:] = N[i,0,:,:]/tinds.size # finish mean

            # Summer
            tinds =  find((Year==years) * ((months==7) + (months==8))) # loop through season
            for tind in tinds:
                salt = d.variables['salt'][tind,:,:,:]
                temp = d.variables['temp'][tind,:,:,:]
                rho = gsw.rho(salt, temp, 0)
                zeta = d.variables['zeta'][tind,:,:]
                zwt = octant.depths.get_zw(d.variables['Vtransform'][:][0], d.variables['Vstretching'][:][0], 
                            salt.shape[0], d.variables['theta_s'][:][0], d.variables['theta_b'][:][0], 
                                h, d.variables['hc'][:][0], zeta=zeta, Hscale=3)
                N[i,1,:,:] = N[i,1,:,:] + np.ma.median(np.sqrt(-g/rho0 * ((rho[2:,:,:]-rho[:-2,:,:])/(zwt[2:,:,:]-zwt[:-2,:,:]))), axis=0) # just save median
            N[i,1,:,:] = N[i,1,:,:]/tinds.size # finish mean

        np.savez(fname2, N=N, xr=grid['xr'].T, yr=grid['yr'].T)

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


    # Horizontal colorbar below plot
    cax = fig.add_axes([0.25, 0.075, 0.5, 0.02]) #colorbar axes
    cb = plt.colorbar(mappable, cax=cax, orientation='horizontal')
    cb.set_label('Mean slope Burger number')

    fig.savefig('figures/slope-burger/' + whichtime + '_mean.png', bbox_inches='tight')


elif whichtime == 'seasonal':

    N = np.zeros((2,f.shape[0],f.shape[1])) # season x grid
 
    if not os.path.exists(fname2):
        # Winter
        i = 0
        tinds =  find(((years>=2004) * (years<=2011)) * ((months==1) + (months==2))) # loop through season
        for tind in tinds:
            salt = d.variables['salt'][tind,:,:,:]
            temp = d.variables['temp'][tind,:,:,:]
            rho = gsw.rho(salt, temp, 0)
            zeta = d.variables['zeta'][tind,:,:]
            zwt = octant.depths.get_zw(d.variables['Vtransform'][:][0], d.variables['Vstretching'][:][0], 
                        salt.shape[0], d.variables['theta_s'][:][0], d.variables['theta_b'][:][0], 
                            h, d.variables['hc'][:][0], zeta=zeta, Hscale=3)
            N[i,:,:] = N[i,:,:] + np.ma.median(np.sqrt(-g/rho0 * ((rho[2:,:,:]-rho[:-2,:,:])/(zwt[2:,:,:]-zwt[:-2,:,:]))), axis=0) # just save median
        N[i,:,:] = N[i,:,:]/tinds.size # finish mean

        # Summer
        i = 1
        tinds =  find((years>=2004 * years<=2011) * ((months==7) + (months==8))) # loop through season
        for tind in tinds:
            salt = d.variables['salt'][tind,:,:,:]
            temp = d.variables['temp'][tind,:,:,:]
            rho = gsw.rho(salt, temp, 0)
            zeta = d.variables['zeta'][tind,:,:]
            zwt = octant.depths.get_zw(d.variables['Vtransform'][:][0], d.variables['Vstretching'][:][0], 
                        salt.shape[0], d.variables['theta_s'][:][0], d.variables['theta_b'][:][0], 
                            h, d.variables['hc'][:][0], zeta=zeta, Hscale=3)
            N[i,:,:] = N[i,:,:] + np.ma.median(np.sqrt(-g/rho0 * ((rho[2:,:,:]-rho[:-2,:,:])/(zwt[2:,:,:]-zwt[:-2,:,:]))), axis=0) # just save median
        N[i,:,:] = N[i,:,:]/tinds.size # finish mean

    np.savez(fname2, N=N, xr=grid['xr'].T, yr=grid['yr'].T)

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


    # Horizontal colorbar below plot
    cax = fig.add_axes([0.25, 0.075, 0.5, 0.02]) #colorbar axes
    cb = plt.colorbar(mappable, cax=cax, orientation='horizontal')
    cb.set_label('Mean slope Burger number')

    fig.savefig('figures/slope-burger/' + whichtime + '_mean.png', bbox_inches='tight')
