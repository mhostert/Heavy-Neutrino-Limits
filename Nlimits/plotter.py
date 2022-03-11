import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
from matplotlib.pyplot import *
from matplotlib.legend_handler import HandlerLine2D
import matplotlib.tri as tri

from . import umu4
from .constraint_dict import *

fsize=11
rcparams={'axes.labelsize':fsize,'xtick.labelsize':fsize,'ytick.labelsize':fsize,\
                'figure.figsize':(1.2*3.7,1.3*2.3617)   }
rc('text', usetex=True)
rc('font',**{'family':'serif', 'serif': ['Computer Modern Roman']})
rcParams.update(rcparams)
axes_form  =[0.15,0.16,0.82,0.76]

label={"e" : r"$|U_{e 4}|^2$", 
       "mu" : r"$|U_{\mu 4}|^2$", 
       "tau" : r"$|U_{\tau 4}|^2$"}

def plot_bound(ax,bound, label='', color='black', lw=0.5, units=1, rasterized=False):
    MN,usqr_bound = bound
    ##############################################
    # Constraints on U\alpha4^2
    l1, = ax.plot(MN*units, usqr_bound, label=label, color=color, lw=lw, rasterized=rasterized)
    ax.fill_between(MN*units, usqr_bound, np.ones(np.size(MN)), 
                    fc='lightgrey', ec='None', lw =0.0, alpha=1,rasterized=rasterized)
    return l1

def plot_ID_individual(axis,listoffiles):
    for i in range(np.size(listoffiles)):
        filename = listoffiles[i]

        if (np.size(filename) == 2):
            name = PATH+filename[0]
        else:
            name = PATH+filename

        m4, Umu4sq = np.genfromtxt(name, unpack=True)
        
        ax.plot(m4/3.0,ualpha4SQR(m4/3.0,gX,m4),label=labels[i],color=colors[i])

def plot_lines(flavor, save=False, invisible=False, m4min=1e-3,m4max=1e2, units = 1):
    fsize = 11
    rc('text', usetex=True)
    params={'axes.labelsize':fsize,'xtick.labelsize':fsize,'ytick.labelsize':fsize,\
                    'figure.figsize':(1.2*3.375,1.5*2.375)  }
    rc('font',**{'family':'serif', 'serif': ['computer modern roman']})
    rcParams.update(params)
    axes_form  = [0.1,0.16,0.85,0.76]


    fig = plt.figure()
    ax = fig.add_axes(axes_form)


    #############
    # get bounds
    list_of_bounds=umu4.get_individual_bounds(muon_bounds, m4min=m4min, m4max=m4max)

    for key in list_of_bounds.keys():
        MN,usqr_bound = list_of_bounds[key]
        ##############################################
        # Constraints on U\alpha4^2
        ax.plot(MN*units, usqr_bound, label=fr"{key.split('_')[0]}")
        ax.fill_between(MN*units, usqr_bound, np.ones(np.size(MN)), 
                        fc='lightgrey', ec='None', lw =0.0, alpha=0.5)



    ax.set_xlim(m4min*units,m4max*units)
    ax.set_ylim(1e-8,2e-2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel(label[flavor])
    ax.set_xlabel(r"$m_{N}$/MeV")
    ax.set_yticks([ 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3,1e-2])
    ax.grid(axis='y', which='both',dashes=(6,1),alpha=0.5,c='black',lw=0.2)
    if save:
        fig.savefig('plots/u'+flavor+'4.pdf')
    return ax, fig