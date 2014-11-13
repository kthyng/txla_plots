'''
Calculate mean wind and variance for wind stress applied in TXLA.
For shelf transport paper.
'''


import numpy as np
import octant
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import netCDF4 as netCDF
from datetime import datetime
from matplotlib.mlab import find


url = 'http://barataria.tamu.edu:6060/thredds/dodsC/NcML/txla_nesting6.nc'

# Times to include in average
year = 2006; 
# startmonth = 7; endmonth = 9; season = 'summer'
startmonth = 1; endmonth = 3; season = 'winter'
startdate = datetime(year, startmonth, 1, 0, 0)
enddate = datetime(year, endmonth, 1, 0, 0)


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
anglev = nc.variables['angle'][:]

x_rho, y_rho = basemap(lon_rho, lat_rho)

u = nc.variables['sustr'][istart:iend, :, :]
v = nc.variables['svstr'][istart:iend, :, :]

# # average
# u = u.mean(axis=0)
# v = v.mean(axis=0)

u = shrink(u, (u.shape[0],mask[1:-1, 1:-1].shape[0],mask[1:-1, 1:-1].shape[1]))
v = shrink(v, (u.shape[0],mask[1:-1, 1:-1].shape[0],mask[1:-1, 1:-1].shape[1]))

u, v = rot2d(u, v, anglev[1:-1, 1:-1])

angle = np.arctan(v/u)
# angle.set_fill_value(np.nan)

np.savez('calcs/wind_stress/calcs' + str(year) + season + '.npz', 
            angle_mean=np.ma.mean(angle), angle_std=np.ma.std(angle))