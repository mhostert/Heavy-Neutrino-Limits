import numpy as np
import numpy.ma as ma
import scipy 

import matplotlib.pyplot as plt
from matplotlib.pyplot import *

from matplotlib import rc, rcParams
import matplotlib.patches as patches
from matplotlib.legend_handler import HandlerLine2D
import matplotlib.tri as tri
from matplotlib import cm
from matplotlib.font_manager import *


fsize=11
rcparams={'axes.labelsize':fsize,'xtick.labelsize':fsize,'ytick.labelsize':fsize,\
                'figure.figsize':(1.2*3.7,1.3*2.3617)   }
rc('text', usetex=True)
rc('font',**{'family':'serif', 'serif': ['Computer Modern Roman']})
rcParams.update(rcparams)
matplotlib.rcParams['hatch.linewidth'] = 0.3

axes_form  =[0.16,0.16,0.81,0.76]
def std_fig(ax_form=axes_form, rasterized=False, figsize=(1.2*3.7,1.3*2.3617)):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes(ax_form, rasterized=rasterized)
    ax.patch.set_alpha(0.0)
    return fig,ax

################################################
def interp_grid(x,y,z, fine_gridx=False, fine_gridy=False, logx=False, logy=False, method='interpolate', smear_stddev=False):

    # default
    if not fine_gridx:
        fine_gridx = 100
    if not fine_gridy:
        fine_gridy = 100

    # log scale x
    if logx:
        xi = np.logspace(np.min(np.log10(x)), np.max(np.log10(x)), fine_gridx)
    else: 
        xi = np.linspace(np.min(x), np.max(x), fine_gridx)
    
    # log scale y
    if logy:
        y = -np.log(y)
        yi = np.logspace(np.min(np.log10(y)), np.max(np.log10(y)), fine_gridy)

    else:
        yi = np.linspace(np.min(y), np.max(y), fine_gridy)

    
    Xi, Yi = np.meshgrid(xi, yi)
    if logy:
        Yi = np.exp(-Yi)

    # triangulation
    if method=='triangulation':
        triang = tri.Triangulation(x, y)
        interpolator = tri.LinearTriInterpolator(triang, z)
        Zi = interpolator(Xi, Yi)
    
    elif method=='interpolate':
        Zi = scipy.interpolate.griddata((x, y), z,\
                                        (xi[None,:], yi[:,None]),\
                                        method='linear', rescale =True)        
    else:
        print(f"Method {method} not implemented.")
    
    # gaussian smear -- not recommended
    if smear_stddev:
            Zi = scipy.ndimage.filters.gaussian_filter(Zi, smear_stddev, mode='nearest', order = 0, cval=0)
    
    return Xi, Yi, Zi


