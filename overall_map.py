'''
Map overall map of various Gulf domains
'''

import numpy as np
from mpl_toolkits.basemap import Basemap, pyproj
import tracpy
import matplotlib.pyplot as plt
import netCDF4 as netCDF
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from scipy.spatial import ConvexHull
import matplotlib

fname = 'TXLA_domain'
doshelf = True
dobay = False
dogulf = False

llcrnrlon=-(108); llcrnrlat=16; 
urcrnrlon=-(74); urcrnrlat=44; projection='lcc'
lat_0=25.; lon_0=-89; resolution='i'; area_thresh=0.
basemap_gulf = Basemap(llcrnrlon=llcrnrlon,
             llcrnrlat=llcrnrlat,
             urcrnrlon=urcrnrlon,
             urcrnrlat=urcrnrlat,
             projection=projection,
             lat_0=lat_0,
             lon_0=lon_0,
             resolution=resolution,
             area_thresh=area_thresh)

if doshelf:
    llcrnrlon=-(98+3/60.+40/3600.); llcrnrlat=25; 
    urcrnrlon=-88; urcrnrlat=30+38/60.+51/3600.; projection='lcc'

    g_shelf = netCDF.Dataset('http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc')

if dobay:
    p = pyproj.Proj(proj='utm', zone='15')

    # Galveston grid
    from sunpy import Grid # suntans code
    p_utm = np.loadtxt('projects/suntans/points_utm.dat')
    lonp, latp = p(p_utm[:,0], p_utm[:,1], inverse=True)
    xp, yp = basemap_gulf(lonp, latp)
    p_lcc = np.zeros(p_utm.shape)
    p_lcc[:,0] = xp
    p_lcc[:,1] = yp
    np.savetxt('projects/suntans/points.dat', p_lcc)
    c_utm = np.loadtxt('projects/suntans/cells_utm.dat')
    c_lcc = c_utm.copy()
    lonp, latp = p(c_utm[:,0], c_utm[:,1], inverse=True)
    xp, yp = basemap_gulf(lonp, latp)
    c_lcc[:,0] = xp
    c_lcc[:,1] = yp
    np.savetxt('projects/suntans/cells.dat', c_lcc)
    grd = Grid('projects/suntans')

if doshelf:
    xr_shelf, yr_shelf = basemap_gulf(g_shelf.variables['lon_rho'][:], g_shelf.variables['lat_rho'][:])
    maskr = g_shelf.variables['mask_rho'][:]

    xr2_shelf = np.ma.masked_where(maskr==0,xr_shelf)
    yr2_shelf = np.ma.masked_where(maskr==0,yr_shelf)

if dogulf:
    g_gulf = netCDF.Dataset('http://barataria.tamu.edu:8080/thredds/dodsC/fmrc/roms/out/ROMS_Output_Feature_Collection_Aggregation_best.ncd')
    xr_gulf2 = np.ma.masked_where(maskr_gulf==0,xr_gulf)
    yr_gulf2 = np.ma.masked_where(maskr_gulf==0,yr_gulf)


### Gulf
fig = plt.figure(figsize=(22,10))
matplotlib.rcParams.update({'font.size': 18})#,'font.weight': 'bold'})
ax_gulf = fig.add_axes([0.001, 0.05, 1.0, 0.9])
basemap_gulf.drawcoastlines(ax=ax_gulf)
basemap_gulf.drawparallels(np.arange(16, 44, 2), dashes=(1, 1), 
                        linewidth=0.15, labels=[1, 0, 0, 0], ax=ax_gulf)
basemap_gulf.drawmeridians(np.arange(-108, -74, 3), dashes=(1, 1), 
                        linewidth=0.15, labels=[0, 0, 0, 1], ax=ax_gulf)
basemap_gulf.drawstates()
basemap_gulf.drawcountries()
basemap_gulf.bluemarble()

if dogulf:
    # Outline numerical domain
    xr_gulf, yr_gulf = basemap_gulf(g_gulf.variables['lon_rho'][:], g_gulf.variables['lat_rho'][:])
    maskr_gulf = g_gulf.variables['mask_rho'][:]

    # Gulf grid
    dx = 1; dy = 1;
    ax_gulf.plot(xr_gulf2[::dx], yr_gulf2[::dy], 
                    xr_gulf2[::dx].T, yr_gulf2[::dy].T, 'grey',alpha=.2)

if doshelf:
    # Shelf grid
    dx = 1; dy = 1;
    ax_gulf.plot(xr2_shelf[::dx], yr2_shelf[::dy], \
                xr2_shelf[::dx].T, yr2_shelf[::dy].T, color='darkcyan',alpha=.1)

if dobay:
    # Plot Galveston grid
    grd.plotmesh(ax=ax_gulf, edgecolors=('grey',), facecolors=('None',), zorder=9)

if dogulf:
    # "Gulf"
    plt.text(0.5, 0.07, 'Gulf', transform = ax_gulf.transAxes, 
            axes=ax_gulf, fontsize=16, color='k', alpha=0.8)


if doshelf:
    ### Shelf

    ax_shelf = zoomed_inset_axes(ax_gulf, 1.6, loc=2) # zoom = 6
    basemap_gulf.drawcoastlines(ax=ax_shelf)
    basemap_gulf.fillcontinents('0.8',ax=ax_shelf)
    basemap_gulf.drawstates()
    basemap_gulf.drawcountries()

if dogulf:
    # Gulf grid
    dx = 1; dy = 1;
    ax_shelf.plot(xr_gulf2[::dx], yr_gulf2[::dy], 
                    xr_gulf2[::dx].T, yr_gulf2[::dy].T, 'lightgrey', alpha=.2, linewidth=.5)

    # Bathymetry
    ax_shelf.contour(xr_gulf, yr_gulf, g_gulf.variables['h'][:], 
                            np.hstack(([10,20],np.arange(50,500,50),np.arange(500,5000,500))), 
                            colors='lightgrey', linewidths=0.5)
if doshelf:
    # Shelf grid
    dx = 1; dy = 1;
    ax_shelf.plot(xr2_shelf[::dx], yr2_shelf[::dy], xr2_shelf[::dx].T, yr2_shelf[::dy].T, 'darkcyan',alpha=.2)

    # sub region of the original image
    x1, x2, y1, y2 = 1080000, 2216090, 630000, 1530000
    ax_shelf.set_xlim(x1, x2)
    ax_shelf.set_ylim(y1, y2)
    plt.xticks(visible=False)
    plt.yticks(visible=False)
    # draw a bbox of the region of the inset axes in the parent axes and
    # connecting lines between the bbox and the inset axes area
    mark_inset(ax_gulf, ax_shelf, loc1=1, loc2=3, fc="none", ec="grey", 
                linewidth=2, zorder=8)

if dobay:
    # Plot Galveston Domain
    # # Get domain polygon from points: http://docs.scipy.org/doc/scipy-dev/reference/generated/scipy.spatial.ConvexHull.html
    # points = np.array([grd.xv,grd.yv]).T
    # hull = ConvexHull(points)
    # for simplex in hull.simplices:
    #   ax_shelf.plot(points[simplex,0], points[simplex,1], color='darkgrey')
    grd.plotmesh(ax=ax_shelf, edgecolors=('grey',), facecolors=('None',), zorder=0)

if doshelf and dobay and dogulf:
    # "Shelf"
    plt.text(0.35, 0.05, 'Shelf', transform = ax_shelf.transAxes, 
            axes=ax_shelf, fontsize=16, color='darkcyan')



if dobay:
    ### Galveston
    ax_bay = zoomed_inset_axes(ax_gulf, 11, loc=1) # zoom = 6
    basemap_gulf.drawcoastlines(ax=ax_bay)
    basemap_gulf.fillcontinents('0.8',ax=ax_bay)
    basemap_gulf.drawstates()
    basemap_gulf.drawcountries()

# # Gulf grid
# ax_bay.plot(xr_gulf2[::dy,::dx], yr_gulf2[::dy,::dx], 
#               xr_gulf2[::dy,::dx].T, yr_gulf2[::dy,::dx].T, 
#               color='black', alpha=.4, linewidth=.5)

if doshelf:
    ### Galveston just area not bay grid
    if not dobay:
        ax_bay = zoomed_inset_axes(ax_gulf, 13, loc=1) # zoom = 6
        basemap_gulf.drawcoastlines(ax=ax_bay)
        basemap_gulf.fillcontinents('0.8',ax=ax_bay)
        basemap_gulf.drawstates()
        basemap_gulf.drawcountries()
        # sub region of the original image
        x1, x2, y1, y2 = 1430000, 1525000, 1300000, 1410000
        ax_bay.set_xlim(x1, x2)
        ax_bay.set_ylim(y1, y2)
        plt.xticks(visible=False)
        plt.yticks(visible=False)
        # draw a bbox of the region of the inset axes in the parent axes and
        # connecting lines between the bbox and the inset axes area
        mark_inset(ax_gulf, ax_bay, loc1=2, loc2=4, fc="none", ec="grey", 
                    linewidth=2, zorder=9)

    # Shelf grid in background of zoom in
    dx = 1; dy = 1;
    ax_bay.plot(xr2_shelf[::dy], yr2_shelf[::dx], xr2_shelf[::dy].T, yr2_shelf[::dx].T, 
        color='darkcyan',alpha=.5, linewidth=.5, zorder=5)

    # Little patch of shelf grid
    ax_bay.plot(xr2_shelf[130:145:dy,265:275:dx], yr2_shelf[130:145:dy,265:275:dx], 
                xr2_shelf[130:145:dy,265:275:dx].T, yr2_shelf[130:145:dy,265:275:dx].T, 
                color='darkcyan', linewidth=1, zorder=5)
    # ax_bay.plot(xr2_shelf[120:128:dy,285:292:dx], yr2_shelf[120:128:dy,285:292:dx], 
    #             xr2_shelf[120:128:dy,285:292:dx].T, yr2_shelf[120:128:dy,285:292:dx].T, 
    #             color='darkcyan', linewidth=1, zorder=5)

if dobay:
    # Galveston grid
    grd.plotmesh(ax=ax_bay, edgecolors=('grey',), facecolors=('None',), zorder=9)

if dogulf:
    # Patch of gulf grid
    ax_bay.plot(xr_gulf2[404:412:dy,126:133:dx], yr_gulf2[404:412:dy,126:133:dx], 
                    xr_gulf2[404:412:dy,126:133:dx].T, yr_gulf2[404:412:dy,126:133:dx].T, 
                    color='black', alpha=1, linewidth=1)

if dobay:
    # sub region of the original image
    x1, x2, y1, y2 = 1420000, 1545000, 1280000, 1410000
    ax_bay.set_xlim(x1, x2)
    ax_bay.set_ylim(y1, y2)
    plt.xticks(visible=False)
    plt.yticks(visible=False)
    # draw a bbox of the region of the inset axes in the parent axes and
    # connecting lines between the bbox and the inset axes area
    mark_inset(ax_gulf, ax_bay, loc1=2, loc2=4, fc="none", ec="grey", 
                linewidth=2, zorder=9)

    # "Bay"
    plt.text(0.05, 0.05, 'Bay', transform = ax_bay.transAxes, 
            axes=ax_bay, fontsize=16, color='grey')


plt.savefig('figures/' + fname + '.png', bbox_inches='tight')