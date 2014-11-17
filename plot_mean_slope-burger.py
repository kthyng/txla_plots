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
from matplotlib import ticker, colors
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
whichseason = 'summer' # 'winter' or 'summer'

loc = 'http://barataria.tamu.edu:6060/thredds/dodsC/NcML/txla_nesting6.nc'
grid = tracpy.inout.readgrid(loc, usebasemap=True)

d = netCDF.Dataset(loc)
t = d.variables['ocean_time'][:]
units = d.variables['ocean_time'].units
dates = netCDF.num2date(t, units)
years = np.asarray([dates[i].year for i in xrange(len(dates))])
months = np.asarray([dates[i].month for i in xrange(len(dates))])
days = np.asarray([dates[i].day for i in xrange(len(dates))])

fname1 = 'calcs/slope-burger/f-alpha.npz'
fname2 = 'calcs/slope-burger/' + whichtime + 'N.npz'

# Things that don't change in time
f = d.variables['f'][:] # Coriolis
pm = d.variables['pm'][:]; pn = d.variables['pn'][:]
h = d.variables['h'][:]
h_smoothed = ndimage.filters.gaussian_filter(h, 20)
dhdy, dhdx = np.gradient(h_smoothed, 1/pn, 1/pm)
alpha = np.sqrt((dhdx)**2 + (dhdy)**2) # dimensional bottom slope, across-shelf
np.savez(fname1, f=f, alpha=alpha, xr=grid['xr'].T, yr=grid['yr'].T)

Years = np.arange(2004,2011)

if whichtime == 'interannual':

    if not os.path.exists(fname2):

        N = np.zeros((Years.size,f.shape[0],f.shape[1])) # season x grid

        for i,Year in enumerate(Years): # loop through years

            # Summer
            count = np.zeros(f.shape) # count of not-masked N values
            tinds =  find((Year==years) * ((months==7) + (months==8))) # loop through season
            for tind in tinds:
                fname = 'calcs/slope-burger/' + 'N-' + str(years[tind]) \
                        + str(months[tind]).zfill(2) + str(days[tind]).zfill(2) + '.npz'
                if not os.path.exists(fname):
                    salt = d.variables['salt'][tind,:,:,:]
                    temp = d.variables['temp'][tind,:,:,:]
                    rho = gsw.rho(salt, temp, 0)
                    zeta = d.variables['zeta'][tind,:,:]
                    zwt = octant.depths.get_zw(d.variables['Vtransform'][:][0], d.variables['Vstretching'][:][0], 
                                salt.shape[0], d.variables['theta_s'][:][0], d.variables['theta_b'][:][0], 
                                    h, d.variables['hc'][:][0], zeta=zeta, Hscale=3)
                    Ntemp = np.ma.median(np.sqrt(-g/rho0 * ((rho[2:,:,:]-rho[:-2,:,:])/(zwt[2:,:,:]-zwt[:-2,:,:]))), axis=0) # just save median
                    Ntemp = Ntemp.filled()
                    np.savez(fname, N=Ntemp) # can't save masked arrays
                else:
                    Ntemp = np.load(fname)['N']
                ind = Ntemp[:,:]>500
                Ntemp[ind] = np.nan
                combined = np.dstack((N[i,:,:],Ntemp))
                N[i,:,:] = np.nansum(combined, axis=2)
                count = count + ~np.isnan(Ntemp)

            N[i,:,:] = N[i,:,:]/count # finish mean
            # pdb.set_trace()

        S = N*alpha/f
        # pdb.set_trace()

        np.savez(fname2, N=N, S=S, f=f, alpha=alpha, xr=grid['xr'].T, yr=grid['yr'].T)
    else:
        S = np.load(fname2)['S']

    # levels = np.linspace(0., 0.5, 11)

    Smean = S.mean(axis=0)

    lev_exp = np.linspace(np.log10(.05), 0.6,11)
    # lev_exp = np.linspace(np.log10(.01), np.ceil(np.log10(np.nanmax(S))),50)
    # lev_exp = np.linspace(np.ceil(np.log10(np.nanmin(S))), np.ceil(np.log10(np.nanmax(S))),50)
    levels = np.power(10, lev_exp)

    fig, axarr = plt.subplots(2,4)
    fig.set_size_inches(13.4, 6.6125)
    fig.subplots_adjust(left=0.03, bottom=0.15, right=1.0, top=0.96, wspace=0.03, hspace=0.11)

    for i, ax in enumerate(axarr.flatten()):
       # Titles for subplots
        if i==4:
            tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2))
            ax.set_title(str(2004+i))
        elif i==7:
            ax.set_frame_on(False)
            ax.set_axis_off()
        else:
            tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2), 
                merslabels=[0, 0, 0, 0], parslabels=[0, 0, 0, 0])
            ax.set_title(str(2004+i))

        if i<7:
            xr=grid['xr'].T; yr=grid['yr'].T
            mappable = ax.contourf(xr, yr, S[i,:,:], cmap='Blues', levels=levels, norm=colors.LogNorm())#, extend='max')
            # mappable = ax.contourf(xr, yr, S[i,:,:]-Smean, cmap='Blues', levels=levels, extend='max')


    # Horizontal colorbar below plot
    cax = fig.add_axes([0.25, 0.075, 0.5, 0.02]) #colorbar axes
    cb = plt.colorbar(mappable, cax=cax, orientation='horizontal')
    cb.set_label('Mean slope Burger number')
    cb.set_ticks(levels[::2])
    cb.set_ticklabels(["%1.2f" % level for level in levels[::2]])

    fig.savefig('figures/slope-burger/' + whichtime + '-' + whichseason + '_mean-log.png', bbox_inches='tight')


elif whichtime == 'seasonal':
 
    if not os.path.exists(fname2):

        N = np.zeros((2,f.shape[0],f.shape[1])) # season x grid
  
        # pdb.set_trace()
        # Winter
        i = 0
        count = np.zeros(f.shape) # count of not-masked N values
        tinds =  find(((years>=2004) * (years<=2010)) * ((months==1) + (months==2))) # loop through season
        for tind in tinds:
            fname = 'calcs/slope-burger/' + 'N-' + str(years[tind]) \
                    + str(months[tind]).zfill(2) + str(days[tind]).zfill(2) + '.npz'
            if not os.path.exists(fname):
                salt = d.variables['salt'][tind,:,:,:]
                temp = d.variables['temp'][tind,:,:,:]
                rho = gsw.rho(salt, temp, 0)
                zeta = d.variables['zeta'][tind,:,:]
                zwt = octant.depths.get_zw(d.variables['Vtransform'][:][0], d.variables['Vstretching'][:][0], 
                            salt.shape[0], d.variables['theta_s'][:][0], d.variables['theta_b'][:][0], 
                                h, d.variables['hc'][:][0], zeta=zeta, Hscale=3)
                Ntemp = np.ma.median(np.sqrt(-g/rho0 * ((rho[2:,:,:]-rho[:-2,:,:])/(zwt[2:,:,:]-zwt[:-2,:,:]))), axis=0) # just save median
                Ntemp = Ntemp.filled()
                np.savez(fname, N=Ntemp) # can't save masked arrays
            else:
                Ntemp = np.load(fname)['N']
            ind = Ntemp[:,:]>500
            Ntemp[ind] = np.nan
            combined = np.dstack((N[i,:,:],Ntemp))
            N[i,:,:] = np.nansum(combined, axis=2)
            count = count + ~np.isnan(Ntemp)

        N[i,:,:] = N[i,:,:]/count # finish mean

        # Summer
        i = 1
        count = np.zeros(f.shape) # count of not-masked N values
        tinds =  find(((years>=2004) * (years<=2010)) * ((months==7) + (months==8))) # loop through season
        for tind in tinds:
            fname = 'calcs/slope-burger/' + 'N-' + str(years[tind]) \
                    + str(months[tind]).zfill(2) + str(days[tind]).zfill(2) + '.npz'
            if not os.path.exists(fname):
                salt = d.variables['salt'][tind,:,:,:]
                temp = d.variables['temp'][tind,:,:,:]
                rho = gsw.rho(salt, temp, 0)
                zeta = d.variables['zeta'][tind,:,:]
                zwt = octant.depths.get_zw(d.variables['Vtransform'][:][0], d.variables['Vstretching'][:][0], 
                            salt.shape[0], d.variables['theta_s'][:][0], d.variables['theta_b'][:][0], 
                                h, d.variables['hc'][:][0], zeta=zeta, Hscale=3)
                Ntemp = np.ma.median(np.sqrt(-g/rho0 * ((rho[2:,:,:]-rho[:-2,:,:])/(zwt[2:,:,:]-zwt[:-2,:,:]))), axis=0) # just save median
                Ntemp = Ntemp.filled()
                np.savez(fname, N=Ntemp) # can't save masked arrays
            else:
                Ntemp = np.load(fname)['N']
            ind = Ntemp[:,:]>500
            Ntemp[ind] = np.nan
            combined = np.dstack((N[i,:,:],Ntemp))
            N[i,:,:] = np.nansum(combined, axis=2)
            count = count + ~np.isnan(Ntemp)

        N[i,:,:] = N[i,:,:]/count # finish mean

        S = N*alpha/f
        # pdb.set_trace()

        np.savez(fname2, N=N, S=S, f=f, alpha=alpha, xr=grid['xr'].T, yr=grid['yr'].T)
    else:
        S = np.load(fname2)['S']


    # Plot Winter and Summer

    lev_exp = np.linspace(np.log10(.05), 0.6,11)
    # lev_exp = np.linspace(np.log10(.01), np.ceil(np.log10(np.nanmax(S))),50)
    # lev_exp = np.linspace(np.ceil(np.log10(np.nanmin(S))), np.ceil(np.log10(np.nanmax(S))),50)
    levels = np.power(10, lev_exp)

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
            # ax.contourf(grid['xr'].T, grid['yr'].T, S[i,:,:], cmap='Blues', levels=levels, extend='max')
            mappable = ax.contourf(grid['xr'].T, grid['yr'].T, S[i,:,:], cmap='Blues', levels=levels, norm=colors.LogNorm())#, extend='max')

        elif i==1:

            tracpy.plotting.background(grid=grid, ax=ax, parslabels=[0,0,0,0], mers=np.arange(-100, -80, 2))
            ax.set_title('Summer')
            # ind = zetas<0
            # mappable = ax.contourf(grid['xr'].T, grid['yr'].T, zetas[ind].mean(axis=0), cmap='Reds', levels=levels, extend='min')
            # mappable = ax.contourf(grid['xr'].T, grid['yr'].T, abs(zetas).mean(axis=0), cmap='Reds', levels=levels, extend='max')
            # mappable = ax.contourf(grid['xr'].T, grid['yr'].T, zetas.mean(axis=0), cmap='RdBu_r', levels=levels, extend='both')
            # mappable = ax.contourf(grid['xr'].T, grid['yr'].T, S[i,:,:], cmap='Blues', levels=levels, extend='max')
            mappable = ax.contourf(grid['xr'].T, grid['yr'].T, S[i,:,:], cmap='Blues', levels=levels, norm=colors.LogNorm())#, extend='max')

        # zmax = np.nanmax((zmaxw,zmaxs))


    # Horizontal colorbar below plot
    cax = fig.add_axes([0.25, 0.075, 0.5, 0.02]) #colorbar axes
    cb = plt.colorbar(mappable, cax=cax, orientation='horizontal')
    cb.set_label('Mean slope Burger number')
    cb.set_ticks(levels[::2])
    cb.set_ticklabels(["%1.2f" % level for level in levels[::2]])

    fig.savefig('figures/slope-burger/' + whichtime + '_mean-log.png', bbox_inches='tight')


    # # Plot seasonal difference
    # levels = np.linspace(0,3.0,11)
    
    # fig = plt.figure(figsize=(6.8375, 6.6125))
    # fig.subplots_adjust(left=0.04, bottom=0.15, right=1.0, top=0.96, wspace=0.07, hspace=0.04)
    # ax = fig.add_subplot(111)
    # tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2))
    # ax.set_title('Winter-Summer Slope Burger Number')
    # mappable = ax.contourf(grid['xr'].T, grid['yr'].T, S[1,:,:]-S[0,:,:], cmap='Blues', levels=levels, extend='max')
    # # Horizontal colorbar below plot
    # cax = fig.add_axes([0.25, 0.075, 0.5, 0.02]) #colorbar axes
    # cb = plt.colorbar(mappable, cax=cax, orientation='horizontal')
    # # cb.set_label('')

    # fig.text(0.125, 0.075, 'Summer', color='#cc0027')
    # fig.text(0.760, 0.075, 'Winter', color='#1b72b7')
