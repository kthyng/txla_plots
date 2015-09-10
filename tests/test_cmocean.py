'''
Test some basic functions of the colormaps.
'''

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def test_cmap_import():
    '''Can we import colormaps and make basic plot.

    '''

    import cmocean

    # find all methods in cmocean
    methods = dir(cmocean)

    # loop through all methods in cmocean
    for method in methods:

        # see if method is a colormap
        if type(eval('cmocean.' + method)) == matplotlib.colors.LinearSegmentedColormap:
            x = np.linspace(0, 10)
            X, _ = np.meshgrid(x, x)
            plt.figure()
            plt.pcolor(X, cmap=eval('cmocean.' + method))
