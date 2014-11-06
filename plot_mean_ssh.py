'''
Plot mean sea surface height in the model domain from the model output. 
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
from matplotlib import ticker
from matplotlib.mlab import find


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


whichtime = 'interannual' # 'seasonal' or 'interannual'
whicharea = 'summer' # 'winter' or 'summer'

loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
grid = tracpy.inout.readgrid(loc, usebasemap=True)

d = netCDF.Dataset(loc)
t = d.variables['ocean_time'][:]
units = d.variables['ocean_time'].units
dates = netCDF.num2date(t, units)
years = np.asarray([dates[i].year for i in xrange(len(dates))])
months = np.asarray([dates[i].month for i in xrange(len(dates))])

fig, axarr = plt.subplots(1,2)
fig.set_size_inches(13.675, 6.6125)
fig.subplots_adjust(left=0.04, bottom=0.15, right=1.0, top=0.96, wspace=0.07, hspace=0.04)

for i, ax in enumerate(axarr):

    # Titles for subplots
    if i==0:
        tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2))
        ax.set_title('Winter')
        # desired time indices for average, 2004-2011 and Jan or Feb
        tinds = (years>=2004) * (years<=2011) * ((months==1) + (months==2))

    elif i==1:
        tracpy.plotting.background(grid=grid, ax=ax, parslabels=[0,0,0,0], mers=np.arange(-100, -80, 2))
        ax.set_title('Summer')
        tinds = (years>=2004) * (years<=2011) * ((months==7) + (months==8))





# fig.savefig('figures/drifters/transport-areas/' + whichtime + '-' + whicharea + '-area.png', bbox_inches='tight')


