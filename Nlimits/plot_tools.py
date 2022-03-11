import scipy 
from scipy.interpolate import splprep, splev

import numpy as np

import colorsys

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
from matplotlib.pyplot import *
from matplotlib.pyplot import cm
import matplotlib.tri as tri
import matplotlib.patches as patches
import matplotlib.colors as mc


def log_interp1d(xx, yy, kind='linear', **kwargs):
    logx = np.log10(xx)
    logy = np.log10(yy)
    lin_interp = scipy.interpolate.interp1d(logx, logy, kind=kind, **kwargs)
    log_interp = lambda zz: np.power(10.0, lin_interp(np.log10(zz)))
    return log_interp


###########################
# Matheus 
fsize=12
fsize_annotate=10

std_figsize = (1.2*3.7,1.3*2.3617)
std_axes_form  =[0.16,0.16,0.81,0.76]

rcparams={'axes.labelsize':fsize,'xtick.labelsize':fsize,'ytick.labelsize':fsize,\
                'figure.figsize':std_figsize, 
                'legend.frameon': False,
                'legend.loc': 'best'  }
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}\usepackage{amssymb}'
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
rcParams.update(rcparams)
matplotlib.rcParams['hatch.linewidth'] = 0.3

# standard figure creation 
def std_fig(ax_form=std_axes_form, 
            figsize=std_figsize,
            rasterized=False):

    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes(ax_form, rasterized=rasterized)
    ax.patch.set_alpha(0.0)

    return fig,ax

# standard saving function
def std_savefig(fig, path, dpi=400, **kwargs):
    fig.savefig(path, dpi = dpi, **kwargs)

# https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib
def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

def double_axes_fig(height = 0.5,
                    gap = 0.1,
                    axis_base = [0.14,0.1,0.80,0.18], 
                    figsize=std_figsize,
                    split_y=False,
                    split_x=False,
                    rasterized=False):

    fig = plt.figure(figsize=figsize)

    if split_y and not split_x:
        axis_base = [0.14,0.1,0.80,0.4-gap/2]
        axis_appended = [0.14,0.5+gap/2,0.80,0.4-gap/2]
    
    elif not split_y and split_x:
        axis_appended = [0.14,0.1,0.4-gap/2,0.8]
        axis_base = [0.14+0.4+gap/2, 0.1, 0.4-gap/2, 0.8]        

    else:
        axis_base[-1] = height
        axis_appended = axis_base+np.array([0, height+gap, 0, 1 - 2*height - gap - axis_base[1] - 0.05])
        

    ax1 = fig.add_axes(axis_appended, rasterized=rasterized)
    ax2 = fig.add_axes(axis_base, rasterized=rasterized)
    ax1.patch.set_alpha(0.0)
    ax2.patch.set_alpha(0.0)

    return fig, [ax1, ax2]

def data_plot(ax, X, Y, xerr, yerr, zorder=2, label='data'):
    return ax.errorbar(X, Y, yerr= yerr, xerr = xerr, \
                    marker="o", markeredgewidth=0.5, capsize=1.0,markerfacecolor="black",\
                    markeredgecolor="black",ms=2,  lw = 0.0, elinewidth=0.8,
                    color='black', label=label, zorder=zorder)

def step_plot(ax, x, y, lw=1, color='red', label='signal', where = 'post', dashes=(3,0), zorder=3):
    return ax.step( np.append(x, np.max(x)+x[-1]),
                    np.append(y, 0.0),
                    where=where,
                    lw = lw, 
                    dashes=dashes,
                    color = color, 
                    label = label, zorder=zorder)


def plot_MB_vertical_region(ax, color='dodgerblue', label=r'MiniBooNE $1 \sigma$'):
    ##########
    # MINIBOONE 2018
    matplotlib.rcParams['hatch.linewidth'] = 0.7
    y = [0,1e10]
    NEVENTS=381.2
    ERROR = 85.2
    xleft = (NEVENTS-ERROR)/NEVENTS
    xright = (NEVENTS+ERROR)/NEVENTS
    ax.fill_betweenx(y,[xleft,xleft],[xright,xright], 
                        zorder=3,
                        ec=color, fc='None',
                        hatch='\\\\\\\\\\',
                        lw=0,
                        label=label)

    ax.vlines(1,0,1e10, zorder=3, lw=1, color=color)
    ax.vlines(xleft,0,1e10, zorder=3, lw=0.5, color=color)
    ax.vlines(xright,0,1e10, zorder=3, lw=0.5, color=color)


def plot_closed_region(points, logx=False, logy=False):
    x,y = points
    if logy:
        if (y==0).any():
            raise ValueError("y values cannot contain any zeros in log mode.")
        sy = np.sign(y)
        ssy = ((np.abs(y)<1)*(-1) + (np.abs(y)>1)*(1))
        y  = ssy*np.log(y*sy)
    if logx:
        if (x==0).any():
            raise ValueError("x values cannot contain any zeros in log mode.")
        sx  = np.sign(x)
        ssx = ((x<1)*(-1) + (x>1)*(1))
        x  = ssx*np.log(x*sx)

    points = np.array([x,y]).T

    points_s     = (points - points.mean(0))
    angles       = np.angle((points_s[:,0] + 1j*points_s[:,1]))
    points_sort  = points_s[angles.argsort()]
    points_sort += points.mean(0)

    tck, u = splprep(points_sort.T, u=None, s=0.0, per=0, k=1)
    u_new = np.linspace(u.min(), u.max(), len(points[:,0]))
    x_new, y_new = splev(u_new, tck, der=0)
    
    if logx:
        x_new = sx*np.exp(ssx*x_new) 
    if logy:
        y_new = sy*np.exp(ssy*y_new) 

    return x_new, y_new


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



def my_histogram(ax, df, var, color='black',label=r'new', density=True, ls='-', var_range=(0,1), units = 1, hatch=''):

    out = ax.hist(df[var, '']*units, 
               bins=30, 
               range=var_range,
               weights=df['weight',], 
                facecolor=color,
              ls=ls,
               edgecolor=color,
               histtype='step',
                density=density,
                 lw=1, zorder=10)

    out = ax.hist(df[var, '']*units, 
               label=label,
               bins=30, 
               range=var_range,
               weights=df['weight',], 
                facecolor='None',
                hatch=hatch,
               edgecolor=color,
                density=density,
                 lw=0.0,
                 alpha =1)    


def error_histogram(ax, df, var, color='black', label=r'new', density=True, ls='-', var_range=(0,1), units = 1, hatch='', cumulative=False):

    w = df['weight',]
    x = df[var, '']*units
    # if cumulative:
    #     var_range = (np.min(df[var, '']), np.max(df[var, '']))

    prediction, bin_edges = np.histogram(x,
             range=var_range,
             bins=50,
             weights=w,
            )

    errors2 = np.histogram(x,
                 range=var_range,
                 bins=50,
                 weights=w**2,
                )[0]
    if cumulative=='cum_sum':
        prediction = np.cumsum(prediction)
        errors2 = np.cumsum(errors2)
    elif cumulative=='cum_sum_prior_to':
        prediction = np.cumsum(prediction[::-1])[::-1]
        errors2 = np.cumsum(errors2[::-1])[::-1]

    area = np.sum(prediction)
    prediction /= area
    errors = np.sqrt(errors2)/area

    # plotting
    ax.plot(bin_edges,
             np.append(prediction, [0]),
             ds='steps-post',
             label=label,
             color=color,
             lw=0.8)

    for edge_left, edge_right, pred, err in zip(bin_edges[:-1], bin_edges[1:], prediction, errors):
        ax.add_patch(
            patches.Rectangle(
            (edge_left, pred-err),
            edge_right-edge_left,
            2 * err,
            color=color,
            hatch=hatch,
            fill=False,
            linewidth=0,
            alpha=0.5,
            )
            )


###########################################
def data_plot(ax, X, DATA, xerr=None, yerr=None, zorder=2, color='black', edgecolor='black'):
    ax.errorbar(X, DATA, yerr= yerr, xerr = xerr, \
                marker="o", markeredgewidth=1.0, capsize=1.0,markerfacecolor=color,\
                markeredgecolor=edgecolor,ms=2, color=color, lw = 0.0, elinewidth=0.8, zorder=zorder)






