
import numpy as np
import octant
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import netCDF4 as netCDF
from datetime import datetime
from matplotlib.mlab import find


tidx = -1       # just get the final frame, for now.
subsample = 9  # roughly every third point in each direction (3**2 = 9)
dx = 25; dy = 25;
scale = 1.0
# url = 'http://testbedapps-dev.sura.org/thredds/dodsC/alldata/Shelf_Hypoxia/tamu/roms/tamu_roms.nc'
url = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'

# Times to include in average
year = 2010; 
startmonth = 7; endmonth = 9; season = 'summer'
# startmonth = 1; endmonth = 3; season = 'winter'
startdate = datetime(year, startmonth, 15, 0, 0)
enddate = datetime(year, endmonth, 15, 0, 0)

iso = 31

# title = '2004'

#####################################################################################

nc = netCDF.Dataset(url)

t = nc.variables['ocean_time'][:]
units = nc.variables['ocean_time'].units
starttime = netCDF.date2num(startdate, units)
endtime = netCDF.date2num(enddate, units)
istart = find(t>=starttime)[0]
iend = find(t>=endtime)[0]

basemap = Basemap(llcrnrlon=-98.5,
                  llcrnrlat=22.5,
                  urcrnrlon=-87.5,
                  urcrnrlat=30.95,
                  projection='lcc',
                  lat_0=30.0,
                  lon_0=-90.0,
                  resolution ='i',
                  area_thresh=0.)
# basemap = Basemap(llcrnrlon=-95.1,
#                   llcrnrlat=27.25,
#                   urcrnrlon=-87.5,
#                   urcrnrlat=30.95,
#                   projection='lcc',
#                   lat_0=30.0,
#                   lon_0=-90.0,
#                   resolution ='i',
#                   area_thresh=0.)


def shrink(a,b):
    """Return array shrunk to fit a specified shape by triming or averaging.
    
    a = shrink(array, shape)
    
    array is an numpy ndarray, and shape is a tuple (e.g., from
    array.shape). a is the input array shrunk such that its maximum
    dimensions are given by shape. If shape has more dimensions than
    array, the last dimensions of shape are fit.
    
    as, bs = shrink(a, b)
    
    If the second argument is also an array, both a and b are shrunk to
    the dimensions of each other. The input arrays must have the same
    number of dimensions, and the resulting arrays will have the same
    shape.
    Example
    -------
    
    >>> shrink(rand(10, 10), (5, 9, 18)).shape
    (9, 10)
    >>> map(shape, shrink(rand(10, 10, 10), rand(5, 9, 18)))        
    [(5, 9, 10), (5, 9, 10)]   
       
    """

    if isinstance(b, np.ndarray):
        if not len(a.shape) == len(b.shape):
            raise Exception, \
                  'input arrays must have the same number of dimensions'
        a = shrink(a,b.shape)
        b = shrink(b,a.shape)
        return (a, b)

    if isinstance(b, int):
        b = (b,)

    if len(a.shape) == 1:                # 1D array is a special case
        dim = b[-1]
        while a.shape[0] > dim:          # only shrink a
            if (dim - a.shape[0]) >= 2:  # trim off edges evenly
                a = a[1:-1]
            else:                        # or average adjacent cells
                a = 0.5*(a[1:] + a[:-1])
    else:
        for dim_idx in range(-(len(a.shape)),0):
            dim = b[dim_idx]
            a = a.swapaxes(0,dim_idx)        # put working dim first
            while a.shape[0] > dim:          # only shrink a
                if (a.shape[0] - dim) >= 2:  # trim off edges evenly
                    a = a[1:-1,:]
                if (a.shape[0] - dim) == 1:  # or average adjacent cells
                    a = 0.5*(a[1:,:] + a[:-1,:])
            a = a.swapaxes(0,dim_idx)        # swap working dim back

    return a

def rot2d(x, y, ang):
    '''rotate vectors by geometric angle'''
    xr = x*np.cos(ang) - y*np.sin(ang)
    yr = x*np.sin(ang) + y*np.cos(ang)
    return xr, yr


mask = nc.variables['mask_rho'][:]
lon_rho = nc.variables['lon_rho'][:]
lat_rho = nc.variables['lat_rho'][:]
h = nc.variables['h'][:]
# angle = nc.variables['angle'][:]
# anglev = nc.variables['angle'][:]

x_rho, y_rho = basemap(lon_rho, lat_rho)

salt = nc.variables['salt'][istart:iend, -1, :, :]
# v = nc.variables['svstr'][istart:iend, :, :]

# average
# salt = salt.mean(axis=0)
# v = v.mean(axis=0)

# u = shrink(u, mask[1:-1, 1:-1].shape)
# v = shrink(v, mask[1:-1, 1:-1].shape)

# u, v = rot2d(u, v, anglev[1:-1, 1:-1])

# # some code to plot random points.
# idx, idy = np.where(mask[1:-1, 1:-1] == 1.0)
# idv = np.arange(len(idx))
# np.random.shuffle(idv)
# Nvec = int(len(idx) / subsample)
# idv = idv[:Nvec]
# idx = idx[idv]
# idy = idy[idv]

figure = plt.figure()
ax = figure.add_subplot(111)

basemap.drawcoastlines()
basemap.fillcontinents()
if year == 2008:
    basemap.drawparallels(np.arange(18, 35), dashes=(1, 1), linewidth=0.15, labels=[1, 0, 0, 0], ax=ax)
    basemap.drawmeridians(np.arange(-100, -80), dashes=(1, 1), linewidth=0.15, labels=[0, 0, 0, 1], ax=ax)
else:
    basemap.drawparallels(np.arange(18, 35), dashes=(1, 1), linewidth=0.15, labels=[0, 0, 0, 0], ax=ax)
    basemap.drawmeridians(np.arange(-100, -80), dashes=(1, 1), linewidth=0.15, labels=[0, 0, 0, 0], ax=ax)

ax.contour(x_rho, y_rho, h, np.hstack(([10,20],np.arange(50,500,50))), 
            colors='lightgrey', linewidths=0.5)

for i in xrange(salt.shape[0]):
    ax.contour(x_rho, y_rho, salt[i,:,:], [iso], colors='0.2', linewidths=0.05)

# q = ax.quiver( x_rho[::dy, ::dx], y_rho[::dy, ::dx], 
#             u[::dy,  ::dx], v[::dy, ::dx], color = '0.3',
#             pivot='middle', zorder=1e35, width=0.003)
#             # scale=1.0/scale, pivot='middle', zorder=1e35, width=0.003)

# if year == 2008:
#     plt.quiverkey(q, 0.85, 0.07, 0.1, label=r'0.1 N m$^{2}$', coordinates='axes')

ax.plot(x_rho[0,:], y_rho[0,:], 'k:')
ax.plot(x_rho[-1,:], y_rho[-1,:], 'k:')
ax.plot(x_rho[:,0], y_rho[:,0], 'k:')
ax.plot(x_rho[:,-1], y_rho[:,-1], 'k:')

plt.title(str(year))

plt.show()

plt.savefig('figures/plume_edge/' + str(year) + season + '.png', bbox_inches='tight', dpi=50)

# plt.close()