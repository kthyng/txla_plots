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
from cmPong import cmPong
from matplotlib.mlab import find
import bisect
from matplotlib import delaunay
import op
import cmocean
import matplotlib.patches as patches
import xarray as xr
import octant


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


start = '2006-01'
stop = '2006-02'



def rot2d(x, y, ang):
    '''rotate vectors by geometric angle'''
    xr = x*np.cos(ang) - y*np.sin(ang)
    yr = x*np.sin(ang) + y*np.cos(ang)
    return xr, yr


# Grid info
loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
# grid_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
# grid = tracpy.inout.readgrid(grid_filename, usebasemap=True, llcrnrlat=22.85, llcrnrlon=-97.9, urcrnrlat=30.5)
# # actually using psi grid here despite the name
# xpsi = np.asanyarray(grid['xpsi'].T, order='C')
# ypsi = np.asanyarray(grid['ypsi'].T, order='C')
# xr = np.asanyarray(grid['xr'].T, order='C')
# yr = np.asanyarray(grid['yr'].T, order='C')
ds = xr.open_dataset(loc)

# current arrows
cdx = 7; cdy = 11 # in indices
wdx = 25; wdy = 40 # in indices, wind arrows

# Colormap for model output
levels = (37-np.exp(np.linspace(0,np.log(37.), 10)))[::-1] # log for salinity, 0 to 36
levels[0] = 0
# levels = (37-np.exp(np.linspace(0,np.log(36.), 10)))[::-1]-1 # log for salinity, 0 to 35
cmap = cmPong.salinity(cmocean.cm.salt, levels)
# cmap = cmPong.salinity('YlGnBu_r', levels)
ilevels = [0,1,2,3,4,5,8] # which levels to label
ticks = [int(tick) for tick in levels[ilevels]] # plot ticks
##

dsgrid = xr.open_dataset('/home/kthyng/shelf/grid.nc')

proj = tracpy.tools.make_proj(setup='nwgom', usebasemap=True)
grid = octant.grid.CGrid_geo(dsgrid['lon_vert'].data, dsgrid['lat_vert'].data, proj)



figname = 'figures/means/salt+vectors' + start + 'to' + stop + '.png'

# Set up plot
fig = plt.figure(figsize=(10.1, 8.4), dpi=100)
ax = fig.add_axes([0.06, 0.00, 0.93, 0.97])
ax.set_frame_on(False)  # kind of like it without the box
proj.drawcoastlines(ax=ax)
proj.fillcontinents('0.8', ax=ax)

proj.drawparallels(np.arange(20,40), dashes=(1, 1), linewidth=0.15, labels=[1, 0, 0, 0], ax=ax)
proj.drawmeridians(np.arange(-100, -80), dashes=(1, 1), linewidth=0.15, labels=[0, 0, 0, 1], ax=ax)
ax.contour(grid.x_rho, grid.y_rho, ds['h'], np.hstack(([10, 20], np.arange(50, 500, 50))), colors='0.2', linewidths=0.5, alpha=0.4)

# tracpy.plotting.background(grid=grid, ax=ax, outline=False, mers=np.arange(-97, -87), merslabels=[0, 0, 1, 0], pars=np.arange(23, 32))

# Label isobaths
ax.text(0.85, 0.865, '10 m', transform=ax.transAxes, fontsize=9, color='0.4', rotation=45)
ax.text(0.88, 0.862, '20', transform=ax.transAxes, fontsize=9, color='0.4', rotation=45)
ax.text(0.87, 0.835, '50', transform=ax.transAxes, fontsize=9, color='0.4', rotation=45)
ax.text(0.89, 0.825, '100', transform=ax.transAxes, fontsize=9, color='0.4', rotation=45)
ax.text(0.9, 0.803, '450', transform=ax.transAxes, fontsize=9, color='0.4', rotation=45)

# # Date
# date = datesModel[itmodel].strftime('%Y %b %02d %H:%M')
# # greyfont = plt.matplotlib.font_manager.FontProperties() # grab default font properties
# # greyfont.set_color('')
# ax.text(0.35, 0.425, date, fontsize=18, color='0.2', transform=ax.transAxes, 
#             bbox=dict(facecolor='white', edgecolor='white', boxstyle='round'))

# Plot surface salinity
# Note: skip ghost cells in x and y so that can properly plot grid cell boxes with pcolormesh
salt = ds['salt'].loc[start:stop].isel(s_rho=-1, eta_rho=slice(1, -1), xi_rho=slice(1, -1)).data.mean(axis=0)
# salt = np.squeeze(m.variables['salt'][itmodel,-1,1:-1,1:-1])
mappable = ax.pcolormesh(grid.x_psi, grid.y_psi, salt, cmap=cmap, vmin=0, vmax=36)

# Surface currents over domain, use psi grid for common locations
u = ds['u'].loc[start:stop].isel(s_rho=-1).data.mean(axis=0)
v = ds['v'].loc[start:stop].isel(s_rho=-1).data.mean(axis=0)
u = op.resize(u, 0)
v = op.resize(v, 1)
anglev = ds['angle'].data
u, v = rot2d(u, v, op.resize(op.resize(anglev, 0), 1))
Q = ax.quiver(grid.x_psi[cdy::cdy,cdx::cdx], grid.y_psi[cdy::cdy,cdx::cdx], u[cdy::cdy,cdx::cdx], v[cdy::cdy,cdx::cdx], 
        color='k', alpha=0.4, pivot='middle', scale=40, width=0.001)
# Q = ax.quiver(xpsi[cdy::cdy,cdy::cdy], ypsi[cdy::cdy,cdy::cdy], Uwind[cdy::cdy,cdy::cdy], Vwind[cdy::cdy,cdy::cdy], 
#         color='k', alpha=0.1, scale=400, pivot='middle', headlength=3, headaxislength=2.8)
qk = ax.quiverkey(Q, 0.18, 0.75, 0.5, r'0.5 m$\cdot$s$^{-1}$ current', labelcolor='0.2', fontproperties={'size': '10'})

# Wind over the domain
# Uwind = w.variables['Uwind'][itwind,:,:]
# Vwind = w.variables['Vwind'][itwind,:,:]
# Uwind, Vwind = rot2d(Uwind, Vwind, anglev)
# Q = ax.quiver(xr[wdy/2::wdy,wdx::wdx], yr[wdy/2::wdy,wdx::wdx], Uwind[wdy/2::wdy,wdx::wdx], Vwind[wdy/2::wdy,wdx::wdx], 
#         color='k', alpha=0.3, scale=300, pivot='middle', headlength=3, headaxislength=2.8)
# qk = ax.quiverkey(Q, 0.18, 0.845, 10, r'10 m$\cdot$s$^{-1}$ wind', labelcolor='0.2', fontproperties={'size': '10'})

sustr = ds['sustr'].loc[start:stop].data.mean(axis=0)
svstr = ds['svstr'].loc[start:stop].data.mean(axis=0)
# sustr = w.variables['sustr'][itwind,:,:]
# svstr = w.variables['svstr'][itwind,:,:]
sustr, svstr = rot2d(op.resize(sustr,1)[1:-1,:], op.resize(svstr,0)[:,1:-1], anglev[1:-1, 1:-1])
Q = ax.quiver(grid.x_rho[wdy+1::wdy,wdx+1::wdx], grid.y_rho[wdy+1::wdy,wdx+1::wdx], sustr[wdy::wdy,wdx::wdx], svstr[wdy::wdy,wdx::wdx], 
        color='k', alpha=0.2, scale=1, pivot='middle', headlength=3, headaxislength=2.8)
qk = ax.quiverkey(Q, 0.18, 0.81, 0.1, label=r'0.0001 N m$^{2}$', labelcolor='0.2', fontproperties={'size': '10'})

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
# plt.close(fig)

