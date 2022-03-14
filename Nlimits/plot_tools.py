import scipy 
from scipy.interpolate import splprep, splev

import numpy as np

import seaborn as sns
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
fsize=12
fsize_annotate=10

std_figsize = (1.2*3.7,1.3*2.3617)
std_axes_form  =[0.16,0.16,0.81,0.76]

# standard figure creation 
def std_fig(ax_form=std_axes_form, 
            figsize=std_figsize,
            rasterized=False):

    rcparams={'axes.labelsize':fsize,'xtick.labelsize':fsize,'ytick.labelsize':fsize,\
                    'figure.figsize':std_figsize, 
                    'legend.frameon': False,
                    'legend.loc': 'best'  }
    plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}\usepackage{amssymb}'
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)
    rcParams.update(rcparams)
    matplotlib.rcParams['hatch.linewidth'] = 0.3
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes(ax_form, rasterized=rasterized)
    ax.patch.set_alpha(0.0)

    return fig,ax

# standard saving function
def std_savefig(fig, path, dpi=400, **kwargs):
    fig.savefig(path, dpi = dpi, **kwargs)
    if '.pdf' in path:
        fig.savefig(path.replace('.pdf','.png'), dpi = dpi, **kwargs)
        fig.savefig(path.replace('.pdf','_white.png'), dpi = dpi, facecolor='white', **kwargs)


def std_plot_limits(case, skip_ids=[], xrange=(5, 1e5), yrange=(1e-10,1e-1), title=None, 
        new_labelpos={}, new_color={}, new_dash={}, grid=False, color_only = [], npoints_interp = 100000, suffix=''):

    fig, ax = std_fig(figsize=(8,4), ax_form=[0.1,0.125,0.88,0.81])
    
    x=np.geomspace(1,1e5, int(npoints_interp))

    # sns.reset_orig()  # get default matplotlib styles back
    labelpos_dic={}
    for i, limit in case.limits.iterrows():
        ilabel, ival = np.nanargmin(limit.interp_func(x)), np.nanmin(limit.interp_func(x))
        labelpos_dic[limit.id] = (x[ilabel]*0.9, ival/2.5)
    color_dic = dict(zip(case.limits['id'],sns.color_palette('tab10', n_colors=len(list(case.limits.iterrows()))))) # a list of RGB tuples
    dash_dic = dict(zip(case.limits['id'], (1+len(color_dic.keys()))*[(1,0)]))
    
    labelpos_dic.update(new_labelpos)
    color_dic.update(new_color)
    dash_dic.update(new_dash)

    for i, limit in case.limits.iterrows():
        if limit.id not in skip_ids:
            
            if len(color_only)>0:

                background_grey = lighten_color('lightgrey', 0.1)
                if limit.id in color_only:
                    c = color_dic[limit.id]
                    LW = 1
                    alpha=1
                else:
                    c = 'black'
                    LW=0.5
                    alpha=1
            else:
                background_grey = lighten_color('lightgrey', 0.3)
                c = color_dic[limit.id]
                alpha=1
                LW = 1

            if (limit.year is None):
                dash = (4,2)
            else:
                dash = dash_dic[limit.id]

            label = fr'\noindent {limit.plot_label}'.replace("(",r" \noindent {\tiny \textsc{(").replace(")", r")} }").replace(r'\\', r"\vspace{-2ex} \\ ")
            ax.plot(x, limit.interp_func(x), zorder=3, color=c, dashes=dash, lw=LW)
            ax.plot(x, limit.interp_func_top(x), color=c, dashes=dash,zorder=2, lw=LW)
            if ('bbn' not in limit.id) and (limit.year is not None):
                ax.fill_between(x, limit.interp_func(x), x/x,  zorder=1, facecolor=background_grey, edgecolor='None', alpha=alpha)    
                # ax.fill_between(x, limit.interp_func(x), limit.interp_func_top(x),  zorder=1, facecolor=background_grey, edgecolor='None', alpha=1)    
            t = ax.annotate(label, xy=labelpos_dic[limit.id], xycoords='data', color=c, zorder=4, fontsize=7.5)
            # t.set_bbox(dict(facecolor=background_grey, alpha=0.2, edgecolor='None'))

    ax.set_yscale("log")
    ax.set_xscale("log")


    ax.set_ylabel(fr"{case.latexflavor}")
    ax.set_xlabel(fr"$m_N/$MeV")

    major=np.array([1e-13,1e-12,1e-11, 1e-10,1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3,1e-2,1e-1,1])
    minor=np.array([2,3,4,5,6,7,8,9])
    minor = np.array([m*minor for m in major]).flatten()[:-9]
    ax.set_yticks(major)
    ax.set_yticks(minor, minor=True)

    if grid:
        ax.grid(axis='y', which='both', dashes=(6,1), alpha=0.25, c='black', lw=0.2)
        ax.grid(axis='x', which='both', dashes=(6,1), alpha=0.25, c='black', lw=0.2)

    ax.set_ylim(*yrange)
    ax.set_xlim(*xrange)
    ax.set_title(title)
    
    std_savefig(fig, path = f'plots/U{case.flavor}N{suffix}.pdf')
    
    return fig, ax 


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






