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
import pdb
import glob


# Times to include in average
year = 2013; 
# startmonth = 7; endmonth = 9; season = 'summer'
startmonth = 1; endmonth = 3; season = 'winter'
startdate = datetime(year, startmonth, 15, 0, 0)
enddate = datetime(year, endmonth, 15, 0, 0)

grid_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
# vert_filename='/atch/raid1/zhangxq/Projects/txla_nesting6/ocean_his_0001.nc'

grid = netCDF.Dataset(grid_filename)


# url = 'http://barataria.tamu.edu:6060/thredds/dodsC/NcML/txla_nesting6.nc'
# url = '/home/kthyng/shelf/' + str(year) + '/ocean_his_*.nc'
if year<=2012:
    years = np.arange(2004, 2013)
    url = []
    for Year in years:
        url.extend(np.sort(glob.glob('/home/kthyng/shelf/' + str(Year) + '/ocean_his_????.nc')))
elif (year==2013) or (year==2014):
    years = np.arange(2013,2015)
    url = []
    for Year in years:
        url.extend(np.sort(glob.glob('/home/kthyng/shelf/' + str(Year) + '/ocean_his_*.nc')))


#####################################################################################

nc = netCDF.MFDataset(url)
# nc = netCDF.Dataset(url)

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


mask = grid.variables['mask_rho'][:]
lon_rho = grid.variables['lon_rho'][:]
lat_rho = grid.variables['lat_rho'][:]
anglev = grid.variables['angle'][:]

x_rho, y_rho = basemap(lon_rho, lat_rho)

# Get wind stress from model for the desired time frame
u = nc.variables['sustr'][istart:iend, :, :]
v = nc.variables['svstr'][istart:iend, :, :]

# Put on the same grid
u = shrink(u, (u.shape[0], mask[1:-1, 1:-1].shape[0], mask[1:-1, 1:-1].shape[1]))
v = shrink(v, (v.shape[0], mask[1:-1, 1:-1].shape[0], mask[1:-1, 1:-1].shape[1]))

# only use wind in a certain area in statistics (winter wind transport region)
lon_rho = shrink(lon_rho, mask[1:-1, 1:-1].shape)
lat_rho = shrink(lat_rho, mask[1:-1, 1:-1].shape)
ind1 = (lon_rho<-90) * (lat_rho>27.5)

# rotate to be on Cartesian grid
u, v = rot2d(u, v, anglev[1:-1, 1:-1])

# Take mean of the stresses in the desired area
u = u[:,ind1].mean()
v = v[:,ind1].mean()

# calculate the angle
angle = np.rad2deg(np.arctan2(v, u))

# Want the break to be at zero degrees
if angle<0:
    angle = angle+360.

np.savez('calcs/wind_stress/calcs' + str(year) + season + '.npz', 
            angle=angle)

