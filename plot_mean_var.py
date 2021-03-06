'''
Make plot of mean variable in color with overlaid mean wind and currents,
for input time period.
'''

import matplotlib as mpl
mpl.use("Agg") # set matplotlib to use the backend that does not require a windowing system
import numpy as np
import os
import netCDF4 as netCDF
import pdb
import matplotlib.pyplot as plt
import tracpy.tools
from datetime import datetime, timedelta
import glob
# from cmPong import cmPong
from matplotlib.mlab import find
import bisect
# from matplotlib import delaunay
import cmocean.cm as cmo
import matplotlib.patches as patches
import xarray as xr
import octant
from matplotlib import cm, colors



# mpl.rcParams['text.usetex'] = True
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


# start = '2006-01'
# stop = '2006-02'

start = '07'  # '07'  # '01'
stop = '08'  # '08'  # '02'
years = [2011]  # np.arange(2004, 2015)
name = None  # None  # 'Summer'  # 'Winter'  # None is to give date in name
var = 'salt'  # 'speed'  # 'salt'
plotcurrents = False  # plot surface currents
plotisohaline = 33  # plot isohaline of value given, or False and don't plot.


if var == 'salt':
    vmin = 0
    vmax = 36
elif var == 'speed':
    vmin = 0
    vmax = 0.5

def rot2d(x, y, ang):
    '''rotate vectors by geometric angle'''
    xr = x*np.cos(ang) - y*np.sin(ang)
    yr = x*np.sin(ang) + y*np.cos(ang)
    return xr, yr

def calc_cmap(cmap=cmo.haline, levels=(37-np.exp(np.linspace(0,np.log(36.), 10)))[::-1]-1):
    '''
    Colormap for salinity for river plumes, with bigger chunks of salinity per color
    section at lower salinity than higher.
    Help from http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps

    Kristen Thyng, Feb 2014

    Inputs:
        cmap        Colormap name to use, e.g. 'YlGnBu'
        levels      edges of colors, as in contourf, to stretch
                    colormap. e.g. for salinity
                    levels = (37-exp(linspace(0,log(36.), 10)))[::-1]-1

    Outputs:
        my_cmap     colormap instance
    '''

    N = levels.size

    # Colors on either side of the edges
    rgb0 = cm.get_cmap(cmap)( np.linspace(0.0, 1.0, N) )[:,0:3]

    red = np.vstack((levels/levels.max(),
                    rgb0[:,0],
                    rgb0[:,0])).T
    red = tuple(map(tuple, red))

    green = np.vstack((levels/levels.max(),
                    rgb0[:,1],
                    rgb0[:,1])).T
    green = tuple(map(tuple, green))

    blue = np.vstack((levels/levels.max(),
                    rgb0[:,2],
                    rgb0[:,2])).T
    blue = tuple(map(tuple, blue))

    cdict = {'red':red, 'green':green, 'blue':blue}
    my_cmap = colors.LinearSegmentedColormap('my_colormap', cdict, 256)

    return my_cmap


# Grid info
loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
# grid_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
# grid = tracpy.inout.readgrid(grid_filename, usebasemap=True, llcrnrlat=22.85, llcrnrlon=-97.9, urcrnrlat=30.5)
# # actually using psi grid here despite the name
# xpsi = np.asanyarray(grid['xpsi'].T, order='C')
# ypsi = np.asanyarray(grid['ypsi'].T, order='C')
# xr = np.asanyarray(grid['xr'].T, order='C')
# yr = np.asanyarray(grid['yr'].T, order='C')
ds = xr.open_dataset(loc, decode_cf=False)
ds['temp'].attrs['missing_value'] = ds['temp'].attrs['_FillValue']
key='salt'
ds[key].attrs['missing_value'] = ds[key].attrs['_FillValue']
ds=xr.decode_cf(ds)

# current arrows
cdx = 7; cdy = 11 # in indices
wdx = 25; wdy = 40 # in indices, wind arrows

# Colormap for model output
if var == 'salt':
    levels = (37-np.exp(np.linspace(0,np.log(37.), 10)))[::-1] # log for salinity, 0 to 36
    levels[0] = 0
    # levels = (37-np.exp(np.linspace(0,np.log(36.), 10)))[::-1]-1 # log for salinity, 0 to 35
    cmap = calc_cmap(cmo.haline, levels)
    # cmap = cmPong.salinity(cmo.haline, levels)
    # cmap = cmPong.salinity('YlGnBu_r', levels)
    ilevels = [0,1,2,3,4,5,8] # which levels to label
    ticks = [int(tick) for tick in levels[ilevels]] # plot ticks
elif var == 'speed':
    cmap = cmo.speed
##

# dsgrid = xr.open_dataset('/atch/raid1/zhangxq/Projects/txla_nesting/txla_grd_v4_new.nc')
dsgrid = xr.open_dataset('../grid.nc')

proj = tracpy.tools.make_proj(setup='nwgom', usebasemap=True)
grid = octant.grid.CGrid_geo(dsgrid['lon_vert'].data, dsgrid['lat_vert'].data, proj)


if name is None:
    figname = 'figures/means/vectors' + var + '-' + start + 'to' + stop + str(years[0]) + '.png'
else:
    figname = 'figures/means/vectors-' + var + '-' + name.lower() + '.png'

# Set up plot
fig = plt.figure(figsize=(10.1, 8.75), dpi=100)
ax = fig.add_axes([0.06, 0.00, 0.93, 0.97])
ax.set_frame_on(False)  # kind of like it without the box
proj.drawcoastlines(ax=ax)
proj.fillcontinents('0.8', ax=ax)

proj.drawparallels(np.arange(20,40), dashes=(1, 1), linewidth=0.15, labels=[1, 0, 0, 0], ax=ax)
proj.drawmeridians(np.arange(-100, -80), dashes=(1, 1), linewidth=0.15, labels=[0, 0, 1, 0], ax=ax)
ax.contour(grid.x_rho, grid.y_rho, ds['h'], np.hstack(([10, 20], np.arange(50, 500, 50))), colors='0.2', linewidths=0.8, alpha=0.5)
ax.contour(grid.x_rho, grid.y_rho, ds['h'], [100], colors='0.0', linewidths=1.5, alpha=0.5)

# tracpy.plotting.background(grid=grid, ax=ax, outline=False, mers=np.arange(-97, -87), merslabels=[0, 0, 1, 0], pars=np.arange(23, 32))

# # Label isobaths
# ax.text(0.85, 0.865, '10 m', transform=ax.transAxes, fontsize=9, color='0.25', rotation=45)
# ax.text(0.88, 0.862, '20', transform=ax.transAxes, fontsize=9, color='0.25', rotation=45)
# ax.text(0.87, 0.835, '50', transform=ax.transAxes, fontsize=9, color='0.25', rotation=45)
# ax.text(0.89, 0.825, '100', transform=ax.transAxes, fontsize=9, color='0.25', rotation=45)
# ax.text(0.9, 0.803, '450', transform=ax.transAxes, fontsize=9, color='0.25', rotation=45)

# Date
# date = datesModel[itmodel].strftime('%Y %b %02d %H:%M')
# # greyfont = plt.matplotlib.font_manager.FontProperties() # grab default font properties
# # greyfont.set_color('')
# date = start.split('-')[0]
if '01' in start and name==None:
    date = 'January-February'
elif '01' in start:
    date = 'Winter'
elif '07' in start and name==None:
    date = 'July-August'
elif '07' in start:
    date = 'Summer'
elif '08' in start:
    date = 'August'
ax.text(0.6, 0.95, date, fontsize=18, color='0.2', transform=ax.transAxes,
            bbox=dict(facecolor='0.8', edgecolor='0.8', boxstyle='round'))

calcsname = 'calcs/means/' + var + '-' + date + str(years[0]) + '.npz'
# import pdb; pdb.set_trace()
if os.path.exists(calcsname):
    d = np.load(calcsname)
    salt = d['salt']; u = d['u']; v = d['v']; sustr = d['sustr']; svstr = d['svstr']
else:
    if var == 'salt':
        # Plot surface salinity
        # Note: skip ghost cells in x and y so that can properly plot grid cell boxes with pcolormesh
        salt = []
        for year in years:
            salt.append(ds[var].loc[str(year) + '-' + start:str(year) + '-' + stop].isel(s_rho=-1, eta_rho=slice(1, -1), xi_rho=slice(1, -1)).data.mean(axis=0))
        salt = np.asarray(salt).mean(axis=0)
        # salt = np.squeeze(m.variables['salt'][itmodel,-1,1:-1,1:-1])
    anglev = ds['angle'].data
    # Surface currents over domain, use psi grid for common locations
    if plotcurrents or var == 'speed':
        u = []
        for year in years:
            u.append(ds['u'].loc[str(year) + '-' + start:str(year) + '-' + stop].isel(s_rho=-1).data.mean(axis=0))
        u = np.asarray(u).mean(axis=0)
        v = []
        for year in years:
            v.append(ds['v'].loc[str(year) + '-' + start:str(year) + '-' + stop].isel(s_rho=-1).data.mean(axis=0))
        v = np.asarray(v).mean(axis=0)
        u = tracpy.op.resize(u, 0)
        v = tracpy.op.resize(v, 1)
        u, v = rot2d(u, v, tracpy.op.resize(tracpy.op.resize(anglev, 0), 1))

    if var == 'speed':
        salt = np.sqrt(u**2 + v**2)

    # wind stress
    sustr = []
    for year in years:
        sustr.append(ds['sustr'].loc[str(year) + '-' + start:str(year) + '-' + stop].data.mean(axis=0))
    sustr = np.asarray(sustr).mean(axis=0)
    svstr = []
    for year in years:
        svstr.append(ds['svstr'].loc[str(year) + '-' + start:str(year) + '-' + stop].data.mean(axis=0))
    svstr = np.asarray(svstr).mean(axis=0)
    # sustr = w.variables['sustr'][itwind,:,:]
    # svstr = w.variables['svstr'][itwind,:,:]
    sustr, svstr = rot2d(tracpy.op.resize(sustr,1)[1:-1,:], tracpy.op.resize(svstr,0)[:,1:-1], anglev[1:-1, 1:-1])
    np.savez(calcsname, salt=salt, u=u, v=v, sustr=sustr, svstr=svstr)

mappable = ax.pcolormesh(grid.x_psi, grid.y_psi, salt, cmap=cmap, vmin=vmin, vmax=vmax)
# import pdb; pdb.set_trace()
if var == 'salt' and plotisohaline:
    ax.contour(grid.x_rho[1:-1,1:-1], grid.y_rho[1:-1,1:-1], salt,
               [plotisohaline], colors='indigo', linewidths=4, alpha=0.4)

if plotcurrents:
    Q = ax.quiver(grid.x_psi[cdy::cdy,cdx::cdx], grid.y_psi[cdy::cdy,cdx::cdx], u[cdy::cdy,cdx::cdx], v[cdy::cdy,cdx::cdx],
            color='k', alpha=0.4, pivot='middle', scale=40, width=0.001)
    # Q = ax.quiver(xpsi[cdy::cdy,cdy::cdy], ypsi[cdy::cdy,cdy::cdy], Uwind[cdy::cdy,cdy::cdy], Vwind[cdy::cdy,cdy::cdy],
    #         color='k', alpha=0.1, scale=400, pivot='middle', headlength=3, headaxislength=2.8)
    qk = ax.quiverkey(Q, 0.18, 0.75, 0.5, r'0.5 m$\cdot$s$^{-1}$ current', labelcolor='0.2', fontproperties={'size': '10'})

Q = ax.quiver(grid.x_rho[wdy+1::wdy,wdx+1::wdx], grid.y_rho[wdy+1::wdy,wdx+1::wdx], sustr[wdy::wdy,wdx::wdx], svstr[wdy::wdy,wdx::wdx],
        color='k', alpha=0.2, scale=1.0, pivot='middle', headlength=3, headaxislength=2.8)
qk = ax.quiverkey(Q, 0.18, 0.81, 0.03, label=r'0.03 N m$^{2}$', labelcolor='0.2', fontproperties={'size': '10'})

# Colorbar in upper left corner
cax = fig.add_axes([0.09, 0.9, 0.35, 0.025]) #colorbar axes
cb = fig.colorbar(mappable, cax=cax, orientation='horizontal')
cb.set_label(r'Surface salinity [g$\cdot$kg$^{-1}$]', fontsize=14, color='0.2')
cb.ax.tick_params(labelsize=14, length=2, color='0.2', labelcolor='0.2')
if var == 'salt':
    cb.set_ticks(ticks)
# # box behind to hide lines
# ax.add_patch( patches.Rectangle( (0.005, 0.925), 0.42, 0.0625, transform=ax.transAxes, color='0.8', zorder=3))
# ax.add_patch( patches.Rectangle( (0.1, 0.895), 0.24, 0.029, transform=ax.transAxes, color='0.8', zorder=3))
# change colorbar tick color http://stackoverflow.com/questions/9662995/matplotlib-change-title-and-colorbar-text-and-tick-colors
cbtick = plt.getp(cb.ax.axes, 'yticklabels')
plt.setp(cbtick, color='0.2')
# pdb.set_trace()

plt.savefig(figname, dpi=300)
plt.close(fig)
