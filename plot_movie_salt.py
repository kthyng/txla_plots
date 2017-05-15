'''
Make salinity plots for movies of the full domain.

To be run in Python 3.
'''

import matplotlib as mpl
mpl.use("Agg") # set matplotlib to use the backend that does not require a windowing system
import numpy as np
import os
import xarray as xr
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from glob import glob
from cmPong import cmPong
from matplotlib.mlab import find
import bisect
from matplotlib import delaunay
import cmocean.cm as cmo
import matplotlib.patches as patches
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import cartopy.feature as cfeature
import pandas as pd


# mpl.rcParams['text.usetex'] = True
mpl.rcParams.update({'font.size': 14})
# mpl.rcParams['font.sans-serif'] = 'Arev Sans, Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Helvetica, Avant Garde, sans-serif'
# mpl.rcParams['mathtext.fontset'] = 'custom'
# mpl.rcParams['mathtext.cal'] = 'cursive'
# mpl.rcParams['mathtext.rm'] = 'sans'
# mpl.rcParams['mathtext.tt'] = 'monospace'
# mpl.rcParams['mathtext.it'] = 'sans:italic'
# mpl.rcParams['mathtext.bf'] = 'sans:bold'
# mpl.rcParams['mathtext.sf'] = 'sans'
# mpl.rcParams['mathtext.fallback_to_cm'] = 'True'


year = 2014



def rot2d(x, y, ang):
    '''rotate vectors by geometric angle'''
    xr = x*np.cos(ang) - y*np.sin(ang)
    yr = x*np.sin(ang) + y*np.cos(ang)
    return xr, yr

# arrows decimation
cdx = 7; cdy = 11 # currents, in indices
wdx = 25; wdy = 40 # wind, in indices

hlevs = [10, 20, 50, 100, 150, 200, 250, 300, 350, 400, 450]

# Grid info
loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_hindcast_agg'
m = xr.open_dataset(loc)

# Rename for convenience
lon_psi = m['lon_psi'][:]
lat_psi = m['lat_psi'][:]
lon_rho = m['lon_rho'][:]
lat_rho = m['lat_rho'][:]
anglev = m.variables['angle'][:]  # theta to rotate wind vectors

## Model output ##
# if year <= 2013:
#     currents_filenames = np.sort(glob.glob('/home/kthyng/shelf/' + str(year) + '/ocean_his_????.nc'))
# elif year == 2014:
#     currents_filenames = np.sort(glob.glob('/home/kthyng/shelf/' + str(year) + '/ocean_his_??.nc'))
# m = netCDF.MFDataset(currents_filenames)

# Time period to use
plotdates = m['ocean_time'].sel(ocean_time=str(year))
# units = m.variables['ocean_time'].units
# starttime = netCDF.date2num(datetime(year, 1, 1, 4, 0, 0), units)
# if year==2014:
#     endtime = netCDF.date2num(datetime(year, 9, 30, 20, 0, 0), units)
# else:
#     endtime = netCDF.date2num(datetime(year+1, 1, 1, 4, 0, 0), units)
# dt = m.variables['ocean_time'][1] - m.variables['ocean_time'][0] # 4 hours in seconds
# ts = np.arange(starttime, endtime, dt)
# itshift = find(starttime==m.variables['ocean_time'][:]) # shift to get to the right place in model output
# datesModel = netCDF.num2date(m.variables['ocean_time'][:], units)
#
# plotdates = netCDF.num2date(ts, units)
# if year == 2014:
#     monthdates = [datetime(year, month, 1, 0, 0, 0) for month in np.arange(1,10)]
# else:
#     monthdates = [datetime(year, month, 1, 0, 0, 0) for month in np.arange(1,13)]

# if not os.path.exists('figures/' + str(year)):
#     os.makedirs('figures/' + str(year))

# Colormap for model output
levels = (37-np.exp(np.linspace(0,np.log(37.), 10)))[::-1] # log for salinity, 0 to 36
levels[0] = 0
# levels = (37-np.exp(np.linspace(0,np.log(36.), 10)))[::-1]-1 # log for salinity, 0 to 35
cmap = cmPong.salinity(cmocean.cm.salt, levels)
# cmap = cmPong.salinity('YlGnBu_r', levels)
ilevels = [0,1,2,3,4,5,8] # which levels to label
ticks = [int(tick) for tick in levels[ilevels]] # plot ticks
##

## Wind forcing ##

# # There are multiple file locations
# if year <= 2012:
#     w = netCDF.Dataset('/atch/raid1/zhangxq/Projects/narr_txla/txla_blk_narr_' + str(year) + '.nc')
# elif year == 2013:
#     w = netCDF.Dataset('/rho/raid/home/kthyng/txla/txla_wind_narr_2013.nc')
# elif year == 2014:
#     w = netCDF.Dataset('/rho/raid/home/kthyng/txla/txla_wind_narr_2014.nc')

# # w = netCDF.Dataset('/atch/raid1/zhangxq/Projects/narr_txla/txla_blk_narr_' + str(year) + '.nc')
# # Wind time period to use
# unitsWind = (w.variables['time'].units).replace('/','-')
# datesWind = netCDF.num2date(w.variables['time'][:], unitsWind)
# # datesWind = datesModel
##

## River forcing ##
Files = sorted(glob('/copano/d1/shared/TXLA_ROMS/inputs/rivers/txla2_river_????_AR_newT_SWpass_weekly.nc'))
ds = [xr.open_dataset(File) for File in Files]
# need to drop extra variable from 2016:
ds[-1] = ds[-1].drop('river_flag')
rds = xr.auto_combine(ds)  # all output here
# take 2/3 of total river inflow as mississippi river discharge
r = (np.abs(rds['river_transport']).sum(axis=1)*2.0/3.0).to_pandas()


# r = xr.open_mfdataset()
# r1 = netCDF.Dataset('/rho/raid/home/kthyng/txla/TXLA_river_4dyes_2012.nc') # use for through 2011
# r2 = netCDF.Dataset('/rho/raid/home/kthyng/txla/TXLA_river_4dyes_2012_2014.nc') # use for 2012-2014
# # River timing
# tr1 = r1.variables['river_time']
# tunitsr1 = tr1.units
# # interpolate times for this data file since at the 12 hours mark instead of beginning of the day
# tr1 = op.resize(tr1, 0)
# datesr1 = netCDF.num2date(tr1[:], tunitsr1)
# tr2 = r2.variables['river_time']
# datesr2 = netCDF.num2date(tr2[:], tr2.units)
# # all of river input
# Q1 = np.abs(r1.variables['river_transport'][:]).sum(axis=1)*2.0/3.0
# # interpolate this like for time
# Q1 = op.resize(Q1, 0)
# Q2 = np.abs(r2.variables['river_transport'][:]).sum(axis=1)*2.0/3.0
# # Combine river info into one dataset
# iend1 = find(datesr1<datetime(2012,1,1,0,0,0))[-1] # ending index for file 1
# tRiver = np.concatenate((tr1[:iend1], tr2[:]), axis=0)
# datesRiver = np.concatenate((datesr1[:iend1], datesr2))
# R = np.concatenate((Q1[:iend1], Q2))
# r1.close(); r2.close()
# # start and end indices in time for river discharge
# itstartRiver = bisect.bisect_left(datesRiver, datetime(year, 1, 1, 0, 0, 0))
# if year == 2014:
#     itendRiver = bisect.bisect_left(datesRiver, datetime(year, 9, 30, 20, 0, 0))
# else:
#     itendRiver = bisect.bisect_left(datesRiver, datetime(year+1, 1, 1, 0, 0, 0))
# ticks for months on river discharge
mticks = [bisect.bisect_left(datesRiver, monthdate) for monthdate in np.asarray(monthdates)]
mticknames = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
##

land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                        edgecolor='face',
                                        facecolor=cfeature.COLORS['land'])


# Loop through times that simulation output is available
for plotdate in plotdates:

    # # Set up before plotting
    # itmodel = bisect.bisect_left(datesModel, plotdate) # index for model output at this time
    # itwind = bisect.bisect_left(datesWind, plotdate) # index for wind at this time
    # itriver = bisect.bisect_left(datesRiver, plotdate) # index for river at this time

    figname = 'figures/salt/movies/' + datesModel[itmodel].isoformat()[0:13] + '.png'

    # Don't redo plot
    if os.path.exists(figname):
        continue

    # Set up plot
    fig = plt.figure(figsize=(9.4, 8.4), dpi=100)
    ax = fig.add_axes([0.06, 0.01, 0.93, 0.95], projection=ccrs.Mercator(central_longitude=-85.0))
    ax.set_frame_on(False) # kind of like it without the box
    ax.set_extent([-98, -87.5, 22.8, 30.5], ccrs.PlateCarree())
    gl = ax.gridlines(linewidth=0.2, color='gray', alpha=0.5, linestyle='-', draw_labels=True)
    # the following two make the labels look like lat/lon format
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    # gl.xlocator = mticker.FixedLocator([-105, -95, -85, -75, -65])  # control where the ticks are
    # gl.xlabel_style = {'size': 15, 'color': 'gray'}  # control how the tick labels look
    # gl.ylabel_style = {'color': 'red', 'weight': 'bold'}
    gl.xlabels_bottom = False  # turn off labels where you don't want them
    gl.ylabels_right = False

    # plot isobaths
    ax.contour(lon_rho, lat_rho, m.h, hlevs, colors='0.6', transform=ccrs.PlateCarree(), linewidths=0.1)

    # Label isobaths
    # ax.text(0.9, 0.92, '10 m', transform=ax.transAxes, fontsize=9, color='0.4', rotation=45)
    ax.text(0.85, 0.865, '10 m', transform=ax.transAxes, fontsize=9, color='0.4', rotation=45)
    ax.text(0.88, 0.862, '20', transform=ax.transAxes, fontsize=9, color='0.4', rotation=45)
    ax.text(0.87, 0.835, '50', transform=ax.transAxes, fontsize=9, color='0.4', rotation=45)
    ax.text(0.89, 0.825, '100', transform=ax.transAxes, fontsize=9, color='0.4', rotation=45)
    ax.text(0.9, 0.803, '450', transform=ax.transAxes, fontsize=9, color='0.4', rotation=45)

    # Date
    datestr = pd.to_datetime(plotdate.data).strftime('%Y %b %02d %H:%M')
    # greyfont = plt.matplotlib.font_manager.FontProperties() # grab default font properties
    # greyfont.set_color('')
    ax.text(0.35, 0.425, datestr, fontsize=18, color='0.2', transform=ax.transAxes,
                bbox=dict(facecolor='white', edgecolor='white', boxstyle='round'))
    # ax.text(0.35, 0.425, date, fontsize=18, color='0.2', transform=ax.transAxes,
    #             bbox=dict(facecolor='white', edgecolor='white', boxstyle='round'))

    # PONG
    ax.text(0.6, 0.97, 'pong.tamu.edu', fontsize=12, transform=ax.transAxes, color='0.3')

    # Plot surface salinity
    # Note: skip ghost cells in x and y so that can properly plot grid cell boxes with pcolormesh
    salt = np.squeeze(m.salt.sel(ocean_time=plotdate).isel(s_rho=-1, eta_rho=slice(1,-1), xi_rho=slice(1,-1)))
    mappable = ax.pcolormesh(lon_psi, lat_psi, salt, cmap=cmap, vmin=0, vmax=36, transform=ccrs.PlateCarree())

    # Mississippi river discharge rate
    axr = fig.add_axes([0.35, 0.2, 0.6, .2])
    # axr.set_frame_on(False) # kind of like it without the box
    for axis in ['top','left','right']:
        axr.spines[axis].set_linewidth(0.05)
    axr.spines['bottom'].set_linewidth(0.0)
    # the plot itself
    axr.fill_between(r[str(year)+'-1-1':datestr].index, r[str(year)+'-1-1':datestr], facecolor='0.4', edgecolor='0.4', zorder=2, alpha=0.5)  # plot up to now
    axr.plot(r[str(year)].index, r[str(year)], '-', color='0.4', alpha=0.3)  # plot whole year
    axr.plot(r[str(year)+'-1-1':datestr].index, r[str(year)+'-1-1':datestr], '-', color='0.4')  # plot up to now
    # horizontal grid lines
    axr.plot([r[str(year)].index[0], r[str(year)].index[-1]], [5, 5], '-', color='0.6', lw=0.5, alpha=0.5)
    axr.plot([r[str(year)].index[0], r[str(year)].index[-1]], [10000, 10000], '-', color='0.6', lw=0.5, alpha=0.5)
    axr.plot([r[str(year)].index[0], r[str(year)].index[-1]], [20000, 20000], '-', color='0.6', lw=0.5, alpha=0.5)
    axr.plot([r[str(year)].index[0], r[str(year)].index[-1]], [30000, 30000], '-', color='0.6', lw=0.5, alpha=0.5)
    # this makes sure alignment stays consistent in different years
    axr.autoscale(axis='x', tight=True)
    axr.set_ylim(-1000,45000)
    # labels
    axr.text(r[str(year)].index[-3]]+16.5, 5, '0', fontsize=9, color='0.4', alpha=0.7)
    axr.text(tRiver[mticks[-3]]+16.5, 10000, '10', fontsize=9, color='0.4', alpha=0.7)
    axr.text(tRiver[mticks[-3]]+16.5, 20000, '20', fontsize=9, color='0.4', alpha=0.7)
    axr.text(tRiver[mticks[-3]]+15, 30000, r'$30\times10^3$ m$^3$s$^{-1}$', fontsize=9, color='0.4', alpha=0.7)
    axr.text(tRiver[mticks[-7]]+15, 30000, 'Mississippi discharge', fontsize=10, color='0.2')
    # ticks
    axr.get_yaxis().set_visible(False)
    axr.get_xaxis().set_visible(False)
    # label month ticks
    for i in xrange(len(mticks)):
        axr.text(tRiver[mticks[i]], 2500, mticknames[i], fontsize=9, color='0.2')
    # make background rectangle so lines don't overlap
    axr.add_patch( patches.Rectangle( (0.3, 0.162), 0.7, 0.2, transform=ax.transAxes, color='white', zorder=1))



    # axr.fill_between(tRiver[itstartRiver:itriver+1], R[itstartRiver:itriver+1], alpha=0.5, facecolor='0.4', edgecolor='0.4', zorder=2)
    # axr.plot(tRiver[itstartRiver:itriver], R[itstartRiver:itriver], '-', color='0.4')
    # axr.plot(tRiver[itriver:itendRiver+1], R[itriver:itendRiver+1], '-', color='0.4', alpha=0.3)
    # axr.plot([tRiver[itstartRiver], tRiver[itendRiver]], [5, 5], '-', color='0.6', lw=0.5, alpha=0.5)
    # axr.plot([tRiver[itstartRiver], tRiver[itendRiver]], [10000, 10000], '-', color='0.6', lw=0.5, alpha=0.5)
    # axr.plot([tRiver[itstartRiver], tRiver[itendRiver]], [20000, 20000], '-', color='0.6', lw=0.5, alpha=0.5)
    # axr.plot([tRiver[itstartRiver], tRiver[itendRiver]], [30000, 30000], '-', color='0.6', lw=0.5, alpha=0.5)
    # this makes sure alignment stays consistent in different years
    axr.autoscale(axis='x', tight=True)
    axr.set_ylim(-1000,45000)
    # labels
    axr.text(tRiver[mticks[-3]]+16.5, 5, '0', fontsize=9, color='0.4', alpha=0.7)
    axr.text(tRiver[mticks[-3]]+16.5, 10000, '10', fontsize=9, color='0.4', alpha=0.7)
    axr.text(tRiver[mticks[-3]]+16.5, 20000, '20', fontsize=9, color='0.4', alpha=0.7)
    axr.text(tRiver[mticks[-3]]+15, 30000, r'$30\times10^3$ m$^3$s$^{-1}$', fontsize=9, color='0.4', alpha=0.7)
    axr.text(tRiver[mticks[-7]]+15, 30000, 'Mississippi discharge', fontsize=10, color='0.2')
    # ticks
    # axr.get_xaxis().set_ticklabels([])
    # axr.xaxis.set_ticks_position('bottom')
    # add ticks for each month
    # axr.set_xticks(tRiver[mticks])
    # axr.tick_params('x', width=1.5, color='0.4') # make ticks bigger
    axr.get_yaxis().set_visible(False)
    axr.get_xaxis().set_visible(False)
    # label month ticks
    for i in xrange(len(mticks)):
        axr.text(tRiver[mticks[i]], 2500, mticknames[i], fontsize=9, color='0.2')
    axr.add_patch( patches.Rectangle( (0.3, 0.162), 0.7, 0.2, transform=ax.transAxes, color='white', zorder=1))

    # Surface currents over domain, use psi grid for common locations
    u = op.resize(np.squeeze(m.variables['u'][itmodel,-1,:,:]), 0)
    v = op.resize(np.squeeze(m.variables['v'][itmodel,-1,:,:]), 1)
    u, v = rot2d(u, v, op.resize(op.resize(anglev, 0), 1))
    Q = ax.quiver(xpsi[cdy::cdy,cdx::cdx], ypsi[cdy::cdy,cdx::cdx], u[cdy::cdy,cdx::cdx], v[cdy::cdy,cdx::cdx],
            color='k', alpha=0.4, pivot='middle', scale=40, width=0.001)
    # Q = ax.quiver(xpsi[cdy::cdy,cdy::cdy], ypsi[cdy::cdy,cdy::cdy], Uwind[cdy::cdy,cdy::cdy], Vwind[cdy::cdy,cdy::cdy],
    #         color='k', alpha=0.1, scale=400, pivot='middle', headlength=3, headaxislength=2.8)
    qk = ax.quiverkey(Q, 0.18, 0.795, 0.5, r'0.5 m$\cdot$s$^{-1}$ current', labelcolor='0.2', fontproperties={'size': '10'})

    # Wind over the domain
    Uwind = w.variables['Uwind'][itwind,:,:]
    Vwind = w.variables['Vwind'][itwind,:,:]
    Uwind, Vwind = rot2d(Uwind, Vwind, anglev)
    Q = ax.quiver(xr[wdy/2::wdy,wdx::wdx], yr[wdy/2::wdy,wdx::wdx], Uwind[wdy/2::wdy,wdx::wdx], Vwind[wdy/2::wdy,wdx::wdx],
            color='k', alpha=0.3, scale=300, pivot='middle', headlength=3, headaxislength=2.8)
    qk = ax.quiverkey(Q, 0.18, 0.845, 10, r'10 m$\cdot$s$^{-1}$ wind', labelcolor='0.2', fontproperties={'size': '10'})

        ax.add_feature(land_10m, facecolor='0.9')

    # sustr = w.variables['sustr'][itwind,:,:]
    # svstr = w.variables['svstr'][itwind,:,:]
    # sustr, svstr = rot2d(op.resize(sustr,1)[1:-1,:], op.resize(svstr,0)[:,1:-1], anglev[1:-1, 1:-1])
    # Q = ax.quiver(xr[wdy+1::wdy,wdx+1::wdx], yr[wdy+1::wdy,wdx+1::wdx], sustr[wdy::wdy,wdx::wdx], svstr[wdy::wdy,wdx::wdx],
    #         color='k', alpha=0.1, scale=4, pivot='middle', headlength=3, headaxislength=2.8)
    # qk = ax.quiverkey(Q, 0.18, 0.65, 0.1, label=r'0.1 N m$^{2}$', labelcolor='0.2', fontproperties={'size': '10'})

    # Colorbar in upper left corner
    cax = fig.add_axes([0.09, 0.91, 0.35, 0.025]) #colorbar axes
    cb = fig.colorbar(mappable, cax=cax, orientation='horizontal')
    cb.set_label(r'Surface salinity [g$\cdot$kg$^{-1}$]', fontsize=14, color='0.2')
    cb.ax.tick_params(labelsize=14, length=2, color='0.2', labelcolor='0.2')
    cb.set_ticks(ticks)
    # box behind to hide lines
    ax.add_patch( patches.Rectangle( (0.005, 0.925), 0.42, 0.0625, transform=ax.transAxes, color='0.8', zorder=3))
    ax.add_patch( patches.Rectangle( (0.1, 0.895), 0.24, 0.029, transform=ax.transAxes, color='0.8', zorder=3))
    # change colorbar tick color http://stackoverflow.com/questions/9662995/matplotlib-change-title-and-colorbar-text-and-tick-colors
    cbtick = plt.getp(cb.ax.axes, 'yticklabels')
    plt.setp(cbtick, color='0.2')
    # pdb.set_trace()

    plt.savefig(figname)
    plt.close(fig)

# To make movie: ffmpeg -r 10 -pattern_type glob -i '2008-*.png' -c:v libx264 -pix_fmt yuv420p -crf 15 movie.mp4
