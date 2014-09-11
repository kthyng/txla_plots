'''
Bathymetry plot of TXLA shelf domain.
'''

# import scipy.io
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
# import visclaw.colormaps as colormaps
from skimage import color
from matplotlib import cm
import matplotlib as mpl
# from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
import tracpy
import tracpy.plotting

mpl.rcParams.update({'font.size': 16})
mpl.rcParams['font.sans-serif'] = 'Arev Sans, Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Helvetica, Avant Garde, sans-serif'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.cal'] = 'cursive'
mpl.rcParams['mathtext.rm'] = 'sans'
mpl.rcParams['mathtext.tt'] = 'monospace'
mpl.rcParams['mathtext.it'] = 'sans:italic'
mpl.rcParams['mathtext.bf'] = 'sans:bold'
mpl.rcParams['mathtext.sf'] = 'sans'
mpl.rcParams['mathtext.fallback_to_cm'] = 'True'


### Plot Bathymetry of TXLA shelf ###

# Read in grid info to get bathy
loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
grid = tracpy.inout.readgrid(loc, usebasemap=True)

sea_cmap = plt.get_cmap('Blues_r')

# levels to plot
levs = np.arange(-100, 10, 10)

# Make plot
fig = plt.figure(figsize=(17,15))
ax = fig.add_subplot(111)
tracpy.plotting.background(grid, fig=fig, ax=ax)
mappable = ax.contourf(grid['xr'], grid['yr'], -grid['h'], cmap=sea_cmap, levels=levs, extend='min')
# ax.set_xlim(lonlims)
# ax.set_ylim(latlims)
# ax.set_xlabel('Longitude [degrees]')
# ax.set_ylabel('Latitude [degrees]')
# # Turn off annoying offset, from https://github.com/clawpack/geoclaw/blob/master/src/python/geoclaw/topotools.py#L844
# ax.ticklabel_format(format="plain", useOffset=False)
# plt.xticks(rotation=20)
cb = fig.colorbar(mappable, shrink=0.775)
cb.set_label('Depth [m]')
# plt.tight_layout()

# Save figure
fig.savefig('figures/domain.png', bbox_inches='tight', pad=0.5)


# ### Plot Bathymetry of Puget Sound, Admiralty Inlet, and Admiralty Head ###

# # download bathymetry, which can be found at: http://figshare.com/preview/_preview/1165560 (27.3MB)

# # Read in bathymetry
# mat = scipy.io.loadmat('cascadia_gridded.mat')

# # x and y limits for these plots
# lonlimsPS = [-123.21, -122.15];
# latlimsPS = [47.02, 48.82];
# lonlimsAI = [-122.8, -122.54]
# latlimsAI = [47.9665, 48.227]
# lonlimsAH = [-122.71, -122.65]
# latlimsAH = [48.12, 48.18]

# # Functionality copied from https://github.com/clawpack/geoclaw/blob/master/src/python/geoclaw/topotools.py#L873
# land_cmap = plt.get_cmap('Greens_r')
# sea_cmap = plt.get_cmap('Blues_r')
# cmapPS = colormaps.add_colormaps((land_cmap, sea_cmap), 
#                                 data_limits=[-325,2500],
#                                 data_break=0.0)
# cmapAI = colormaps.add_colormaps((land_cmap, sea_cmap), 
#                                 data_limits=[-200,175],
#                                 data_break=0.0)
# cmapAH = colormaps.add_colormaps((land_cmap, sea_cmap), 
#                                 data_limits=[-110,50],
#                                 data_break=0.0)

# # levels to plot
# levsPS = np.concatenate((np.arange(-325, 0, 25), np.arange(0,3000,500)))
# levsAI = np.concatenate((np.arange(-200, 0, 20), np.arange(0,350,175))) #200,25)))
# levsAH = np.concatenate((np.arange(-120, 0, 20), np.arange(0,100,50)))


# # Make Puget Sound plot
# fig = plt.figure(figsize=(16,16))
# axPS = fig.add_subplot(111)
# mappablePS = axPS.contourf(mat['lon_topo'], mat['lat_topo'], mat['z_topo'], cmap=cmapPS, levels=levsPS)
# # outline coast in case plot is printed
# axPS.contour(mat['lon_topo'], mat['lat_topo'], mat['z_topo'], [0], lw=3, colors='0.15')
# axPS.set_xlim(lonlimsPS)
# axPS.set_ylim(latlimsPS)
# axPS.set_xlabel('Longitude [degrees]')
# axPS.set_ylabel('Latitude [degrees]')
# # Turn off annoying offset, from https://github.com/clawpack/geoclaw/blob/master/src/python/geoclaw/topotools.py#L844
# axPS.ticklabel_format(format="plain", useOffset=False)
# plt.xticks(rotation=20)
# cbPS = fig.colorbar(mappablePS)
# cbPS.set_label('Height/depth [m]')
# plt.tight_layout()
# # Label
# axPS.text(0.7, 0.025, 'Puget Sound', transform=axPS.transAxes, color='0.15')


# # Inset magnified plot of Admiralty Inlet
# axAI = zoomed_inset_axes(axPS, 2, loc=1)
# mappableAI = axAI.contourf(mat['lon_topo'], mat['lat_topo'], mat['z_topo'], cmap=cmapAI, levels=levsAI)
# # outline coast in case plot is printed
# axAI.contour(mat['lon_topo'], mat['lat_topo'], mat['z_topo'], [0], lw=3, colors='0.15')
# axAI.set_xlim(lonlimsAI)
# axAI.set_ylim(latlimsAI)
# # turn off ticks
# plt.xticks(visible=False)
# plt.yticks(visible=False)
# plt.setp(axAI,xticks=[],yticks=[])
# # Inlaid colorbar
# caxAI = fig.add_axes([0.735, 0.77, 0.0125, 0.2])
# cbAI = plt.colorbar(mappableAI, cax=caxAI, orientation='vertical')
# # draw a bbox of the region of the inset axes in the parent axes and
# # connecting lines between the bbox and the inset axes area
# mark_inset(axPS, axAI, loc1=2, loc2=4, fc="none", ec="0.3", lw=1.5)
# # Label
# axAI.text(0.044, 0.04, 'Admiralty Inlet', transform=axAI.transAxes, color='0.15')


# # Inset magnified plot of Admiralty Head
# axAH = zoomed_inset_axes(axPS, 8, loc=3)
# mappableAH = axAH.contourf(mat['lon_topo'], mat['lat_topo'], mat['z_topo'], cmap=cmapAH, levels=levsAH)
# # outline coast in case plot is printed
# axAH.contour(mat['lon_topo'], mat['lat_topo'], mat['z_topo'], [0], lw=3, colors='0.15')
# axAH.set_xlim(lonlimsAH)
# axAH.set_ylim(latlimsAH)
# # turn off ticks
# plt.xticks(visible=False)
# plt.yticks(visible=False)
# plt.setp(axAH,xticks=[],yticks=[])
# # Inlaid colorbar
# caxAH = fig.add_axes([0.35, 0.1, 0.0125, 0.2])
# cbAH = plt.colorbar(mappableAH, cax=caxAH, orientation='vertical')
# # draw a bbox of the region of the inset axes in the parent axes and
# # connecting lines between the bbox and the inset axes area
# mark_inset(axPS, axAH, loc1=2, loc2=4, fc="none", ec="0.3", lw=1.5)
# # Label
# axAH.text(0.45, 0.92, 'Admiralty Head', transform=axAH.transAxes, color='0.15')

# plt.draw()
# plt.show()

# # Save figure
# fig.savefig('figures/domains.png')
# 