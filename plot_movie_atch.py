'''
Make Mississippi river dye plots for movies of the full domain.

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
import cmocean.cm as cmo
import matplotlib.patches as patches
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import cartopy.feature as cfeature
import pandas as pd
from matplotlib import cm, colors
mpl.rcParams.update({'font.size': 12})

def rot2d(x, y, ang):
    '''rotate vectors by geometric angle'''
    xr = x*np.cos(ang) - y*np.sin(ang)
    yr = x*np.sin(ang) + y*np.cos(ang)
    return xr, yr

def resize(A, dim):
    """
    Average neighboring elements in an array A along a dimension dim.

    Args:
        A (array): array size, [m x n] in 2D case. Can be up to 3D.
        dim (int): dimension on which to act. Can be up to 2 (0, 1, 2).

    Returns:
        * A - array of size [(m-1) x n] if dim is 0
    """

    # B is A but rolled such that the dimension that is to be resized is in
    # the 0 position
    B = np.rollaxis(A, dim)

    # Do averaging
    B = 0.5*(B[0:-1]+B[1:])

    # Roll back to original
    return np.rollaxis(B, 0, dim+1)


# arrows decimation
cdx = 7; cdy = 11 # currents, in indices
wdx = 25; wdy = 40 # wind, in indices

hlevs = [10, 20, 50, 100, 150, 200, 250, 300, 350, 400, 450]  # isobath contour depths

# Grid info
try:
    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_hindcast_agg'
    m = xr.open_dataset(loc)
except:
    loc = 'http://copano.tamu.edu:8080/thredds/dodsC/NcML/txla_hindcast_agg'
    m = xr.open_dataset(loc)
#loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_hindcast_agg'
#m = xr.open_dataset(loc)

# Rename for convenience
lon_psi = m['lon_psi'][:].data
lat_psi = m['lat_psi'][:].data
lon_rho = m['lon_rho'][:].data
lat_rho = m['lat_rho'][:].data
anglev = m.variables['angle'][:].data  # theta to rotate wind vectors

# Colormap for model output
cmap = cmo.tempo
cmin = 0.1; cmax = 100; dc = 10
# ticks = [0.1, 1, 10, 100]  # np.arange(cmin, cmax+dc, dc)
# ticklabels = ['0.1', '1', '10', '100']

varname = 'dye_03'
label = r'Atchafalaya river dye [% river water]'
factor = 100  # to change to percentage

os.makedirs('figures/%s/movies' % varname, exist_ok=True)

land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                        edgecolor='face',
                                        facecolor=cfeature.COLORS['land'])
states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='10m',
    facecolor='none')
river = cfeature.NaturalEarthFeature('physical', 'rivers_lake_centerlines', '10m',
                                        edgecolor='b',
                                        facecolor='none')
rivermore = cfeature.NaturalEarthFeature('physical', 'rivers_north_america', '10m')

## River forcing ##
try:
    Files = sorted(glob('/copano/d1/shared/TXLA_ROMS/inputs/rivers/txla2_river_????_AR_newT_SWpass_weekly.nc'))
    ds = [xr.open_dataset(File) for File in Files]
    # need to drop extra variable from 2016:
    ds[-1] = ds[-1].drop('river_flag')
except:
    # in case I am running on rainier with expandrive
    Files = sorted(glob('/Volumes/copano.tamu.edu/d1/shared/TXLA_ROMS/inputs/rivers/txla2_river_????_AR_newT_SWpass_weekly.nc'))
    # Files.pop(-1)  # have to remove 2016 because the file isn't working
    ds = [xr.open_dataset(File) for File in Files]
    # # need to drop extra variable from 2016:
    # ds[-1] = ds[-1].drop('river_flag')
rds = xr.auto_combine(ds)  # all output here
# take 2/3 of total river inflow as mississippi river discharge
r = (np.abs(rds['river_transport']).sum(axis=1)*2.0/3.0).to_pandas()

base = 'figures/' + varname + '/movies/'
years = np.arange(1993, 2018)

for year in years:

    # Time period to use
    plotdates = m['ocean_time'].sel(ocean_time=str(year))

    mticknames = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
    iticks = [datetime(year, month, 2, 0, 0) for month in np.arange(1,13)]


    # Loop through times that simulations were started
    for plotdate in plotdates:

        figname = base + pd.to_datetime(plotdate.data).isoformat()[0:13] + '.png'

        # Don't redo plot
        if os.path.exists(figname):
            continue

        # Set up plot
        fig = plt.figure(figsize=(9.4, 7.7), dpi=100)
        ax = fig.add_axes([0.06, 0.01, 0.93, 0.95], projection=ccrs.Mercator(central_longitude=-85.0))
        ax.set_frame_on(False) # kind of like it without the box
        ax.set_extent([-98, -87.5, 22.8, 30.5], ccrs.PlateCarree())
        gl = ax.gridlines(linewidth=0.2, color='gray', alpha=0.5, linestyle='-', draw_labels=True)
        # the following two make the labels look like lat/lon format
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabels_bottom = False  # turn off labels where you don't want them
        gl.ylabels_right = False

        # plot isobaths
        ax.contour(lon_rho, lat_rho, m.h, hlevs, colors='0.6', transform=ccrs.PlateCarree(), linewidths=0.5)

        # Date
        datestr = pd.to_datetime(plotdate.data).strftime('%Y %b %d %H:%M')
        ax.text(0.35, 0.425, datestr, fontsize=18, color='0.2', transform=ax.transAxes,
                    bbox=dict(facecolor='white', edgecolor='white', boxstyle='round'))

        # PONG
        ax.text(0.55, 0.97, 'pong.tamu.edu', fontsize=12, transform=ax.transAxes, color='0.3')

        # Plot variable
        # Note: skip ghost cells in x and y so that can properly plot grid cell boxes with pcolormesh
        var = m[varname].sel(ocean_time=plotdate).isel(s_rho=-1, eta_rho=slice(1,-1), xi_rho=slice(1,-1))
        mappable = ax.pcolormesh(lon_psi, lat_psi, var*factor, cmap=cmap,
                                 vmin=cmin, vmax=cmax,
                                 transform=ccrs.PlateCarree(),
                                 norm=colors.LogNorm(vmin=cmin, vmax=cmax))
        ax.add_feature(land_10m, facecolor='0.8')
        ax.add_feature(river, edgecolor='cornflowerblue', facecolor='none', linewidth=1)
        ax.add_feature(rivermore, edgecolor='cornflowerblue', facecolor='none', linewidth=0.5)
        ax.coastlines(resolution='10m')  # coastline resolution options are '110m', '50m', '10m'
        ax.add_feature(states_provinces, edgecolor='0.2')
        ax.add_feature(cfeature.BORDERS, linestyle='-', edgecolor='0.2')

        # Mississippi river discharge rate
        axr = fig.add_axes([0.33, 0.2, 0.6, .2])
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
        axr.text(r[str(year)+'-10-15'].index, 5, '0', fontsize=9, color='0.4', alpha=0.7)
        axr.text(r[str(year)+'-10-15'].index, 10000, '10', fontsize=9, color='0.4', alpha=0.7)
        axr.text(r[str(year)+'-10-15'].index, 20000, '20', fontsize=9, color='0.4', alpha=0.7)
        axr.text(r[str(year)+'-10-15'].index, 30000, r'30$\times$10$^3$ m$^3$s$^{-1}$', fontsize=9, color='0.4', alpha=0.7)
        axr.text(r[str(year)+'-6-15'].index, 30000, 'Mississippi discharge', fontsize=9, color='0.4', alpha=0.7)
        # ticks
        axr.get_yaxis().set_visible(False)
        axr.get_xaxis().set_visible(False)
        # label month ticks
        [axr.text(itick, 2500, name, fontsize=9, color='0.2') for itick, name in zip(iticks, mticknames)]

        # Surface currents over domain, use psi grid for common locations
        u = resize(m.u.sel(ocean_time=plotdate).isel(s_rho=-1).data, 0)
        v = resize(m.v.sel(ocean_time=plotdate).isel(s_rho=-1).data, 1)
        u, v = rot2d(u, v, resize(resize(anglev, 0), 1))
        Q = ax.quiver(lon_psi[cdy::cdy,cdx::cdx], lat_psi[cdy::cdy,cdx::cdx], u[cdy::cdy,cdx::cdx], v[cdy::cdy,cdx::cdx],
                color='k', alpha=0.4, pivot='middle', scale=40, width=0.001, transform=ccrs.PlateCarree())
        qk = ax.quiverkey(Q, 0.11, 0.79, 0.5, r'0.5 m$\cdot$s$^{-1}$ current', labelcolor='0.2', fontproperties={'size': '10'})

        # Wind over the domain
        Uwind = m.Uwind.sel(ocean_time=plotdate).data
        Vwind = m.Vwind.sel(ocean_time=plotdate).data
        Uwind, Vwind = rot2d(Uwind, Vwind, anglev)
        wdy2 = int(wdy/2)
        Q = ax.quiver(lon_rho[wdy2::wdy,wdx::wdx], lat_rho[wdy2::wdy,wdx::wdx], Uwind[wdy2::wdy,wdx::wdx], Vwind[wdy2::wdy,wdx::wdx],
                color='k', alpha=0.3, scale=300, pivot='middle', headlength=3, headaxislength=2.8, transform=ccrs.PlateCarree())
        qk = ax.quiverkey(Q, 0.11, 0.84, 10, r'10 m$\cdot$s$^{-1}$ wind', labelcolor='0.2', fontproperties={'size': '10'})

        # Colorbar in upper left corner
        cax = fig.add_axes([0.08, 0.92, 0.32, 0.018]) #colorbar axes
        cb = fig.colorbar(mappable, cax=cax, orientation='horizontal')
        cb.set_label(label, fontsize=13, color='0.2')
        cb.ax.tick_params(labelsize=12, length=2, color='0.2', labelcolor='0.2')
        # cb.set_ticks(ticks)
        # cb.ax.set_yticklabels(ticklabels)  # horizontal colorbar
        # change colorbar tick color http://stackoverflow.com/questions/9662995/matplotlib-change-title-and-colorbar-text-and-tick-colors
        cbtick = plt.getp(cb.ax.axes, 'yticklabels')
        plt.setp(cbtick, color='0.2')
        plt.savefig(figname)
        plt.close(fig)

    if not os.path.exists(base + str(year) + '_low.mp4'):
        # low resolution animation
        command = 'ffmpeg -r 15 -pattern_type glob -i "' + base + str(year) + '-*.png" -c:v libx264 -pix_fmt yuv420p -crf 30 ' + base + str(year) + '_low.mp4'
        os.system(command)
    if not os.path.exists(base + str(year) + '_high.mp4'):
        # high resolution animation
        command = 'ffmpeg -r 15 -pattern_type glob -i "' + base + str(year) + '-*.png" -c:v libx264 -pix_fmt yuv420p -crf 20 ' + base + str(year) + '_high.mp4'
        os.system(command)
