'''
Use this to make surface salinity near river plume have an appropriately
stretched colorbar given the way the surface area grows.

Suggested usage:
import cmstretch
cmap = cmstretch.cmap()
mappable = ax.pcolormesh(x, y, salt, cmap=cmap, vmin=0, vmax=36)
cb = fig.colorbar(mappable)
cb.set_ticks(ticks)

'''

import cmocean.cm as cmo
import numpy as np
from matplotlib import cm, colors

def cmap(cmap=cmo.haline, levels=(37-np.exp(np.linspace(0,np.log(36.), 10)))[::-1]-1):
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
    rgb0 = cm.get_cmap(cmap)(np.linspace(0.0, 1.0, N))[:,0:3]

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

    # return information for labeling ticks nicely since stretched
    ilevels = [0,1,2,3,4,5,8] # which levels to label
    ticks = [int(tick) for tick in levels[ilevels]] # plot ticks

    return my_cmap, ticks

if __name__ == "__main__":

    import netCDF4 as netCDF
    import matplotlib.pyplot as plt

    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/fmrc/txla_10yr_mpdata_his/out/TXLA_10yr_MPDATA_HIS_best.ncd'
    d = netCDF.Dataset(loc)
    salt = d['salt'][0,-1,:,:]
    lonr = d['lon_rho'][:]
    latr = d['lat_rho'][:]
    cmap, ticks = cmap()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    mappable = ax.pcolormesh(lonr, latr, salt, cmap=cmap, vmin=0, vmax=36)
    cb = fig.colorbar(mappable)
    cb.set_ticks(ticks)
    plt.show()
