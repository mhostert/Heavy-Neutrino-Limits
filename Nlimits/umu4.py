import numpy as np
from matplotlib import rc, rcParams
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from scipy import interpolate
import scipy.stats
import sys
import scipy.ndimage as ndimage
import os
from scipy import interpolate

from .constraint_dict import *



## Return a function that takes m4 [GeV] and gives bound on Ualpha4^2
def get_individual_bounds(dict_bounds, invisible=False, m4min=1e-3, m4max=1e4):
    GREYCOLOR = "grey"
    ALP = 0.5
    y=[]
    res_lines = {}

    for key in dict_bounds.keys():
       
        bound = dict_bounds[key]

        ## FIX ME -- no functionality for gaps between top of constraint and other bounds
        # check if it is a closed contour with top and bottom files
        if (np.size(bound['file']) == 2):
            this_file = bound['file'][0]
        else:
            this_file = bound['file']

        m4, Umu4sq = np.genfromtxt(this_file, unpack=True)

        # make sure the data points are ordered in m4
        order = np.argsort(m4)
        m4 = np.array(m4)[order]
        Umu4sq = np.array(Umu4sq)[order]

        if bound['units'] == 'GeV':
            units = 1
        elif bound['units'] == 'MeV':
#             axis.plot(m4, Umu4sq, c=colors[i], ls =linestyles[i], lw=2.0, label=labels[i])
            units = 1e-3
        elif bound['units'] == 'keV':
#             axis.plot(m4, Umu4sq, c=colors[i], ls =linestyles[i], lw=2.0, label=labels[i])
            units = 1e-6
        elif bound['units'] == 'eV':
#             axis.plot(m4, Umu4sq, c=colors[i], ls =linestyles[i], lw=2.0, label=labels[i])
            units = 1e-9
        else:
            print(f"Skipping {key} bound with unknown units of {bound['units']}.")

        if ~(np.min(m4*units) > m4max or np.max(m4*units) < m4min):
            # support of the bound
            x=np.logspace(np.log10(np.min(m4*units)),np.log10(np.max(m4*units)),1000)

            # Only append if invisible
            if invisible:
                if bound['invisible']:
                    f = interpolate.interp1d(m4*units, Umu4sq, kind='linear', bounds_error=False, fill_value=1.0, assume_sorted=False)    
                    res_lines[key] = [x, f(x)]
            # append all bounds
            else: 
                f = interpolate.interp1d(m4*units, Umu4sq, kind='linear', bounds_error=False, fill_value=1.0, assume_sorted=False)    
                res_lines[key] = [x, f(x)]

    return res_lines



## Return a function that takes m4 [GeV] and gives bound on Ualpha4^2
def get_combined_bounds(dict_bounds, invisible=False):
    GREYCOLOR = "grey"
    ALP = 0.5
    y=[]

    x=np.logspace(-3,2,10000)
    for key in dict_bounds.keys():
       
        bound = dict_bounds[key]

        ## FIX ME -- no functionality for gaps between top of constraint and other bounds
        # check if it is a closed contour with top and bottom files
        if (np.size(bound['file']) == 2):
            this_file = bound['file'][0]
        else:
            this_file = bound['file']

        m4, Usqr = np.genfromtxt(this_file, unpack=True)

        # make sure the data points are ordered in m4
        order = np.argsort(m4)
        m4 = np.array(m4)[order]
        Usqr = np.array(Usqr)[order]


        if bound['units'] == 'GeV':
            units = 1
        elif bound['units'] == 'MeV':
            units = 1e-3
        elif bound['units'] == 'keV':
            units = 1e-6
        elif bound['units'] == 'eV':
            units = 1e-9
        else:
            print(f"Skipping {key} bound with unknown units of {bound['units']}.")


        # Only append if invisible
        if invisible:
            if bound['invisible']:
                f = interpolate.interp1d(m4*units, Usqr, kind='linear', bounds_error=False, fill_value=1.0, assume_sorted=False)    
                y.append(f(x))
        # append all bounds
        else: 
            f = interpolate.interp1d(m4*units, Usqr, kind='linear', bounds_error=False, fill_value=1.0, assume_sorted=False)    
            y.append(f(x))

    y = np.array(y)
    z = np.ones(np.size(y[0]))


    for i in range(0,np.shape(y)[0]):
        for j in range(0, np.size(y[i])):
            if y[i,j] < z[j]:
                z[j] = y[i,j]

    return interpolate.interp1d(x, z, kind='linear', bounds_error=False, fill_value=1.0, assume_sorted=False)    


def get_interp(bound, this_file, invisible=False):
    

    if bound['units'] == 'GeV':
        units = 1
    elif bound['units'] == 'MeV':
        units = 1e-3
    elif bound['units'] == 'keV':
        units = 1e-6
    elif bound['units'] == 'eV':
        units = 1e-9
    else:
        print(f"Skipping {key} bound with unknown units of {bound['units']}.")


    m4, Usqr = np.genfromtxt(this_file, unpack=True)

    # # make sure the data points are ordered in m4
    # order = np.argsort(m4)
    # m4 = np.array(m4)[order]
    # Usqr = np.array(Usqr)[order]

    # Only append if invisible
    if invisible:
        if bound['invisible']:
            return interpolate.interp1d(m4*units, Usqr, kind='linear', bounds_error=False, fill_value=np.nan, assume_sorted=False)    
        # append all bounds
        else: 
            return interpolate.interp1d(m4*units, Usqr, kind='linear', bounds_error=False, fill_value=np.nan, assume_sorted=False)    



## Return a function that takes m4 [GeV] and gives bound on Ualpha4^2
def plot_individual(ax, dict_bounds, invisible=False):
    GREYCOLOR = "grey"
    ALP = 0.5
    y=[]

    x=np.logspace(-3,3,1000)
    for key in dict_bounds.keys():
       
        bound = dict_bounds[key]

        ## FIX ME -- no functionality for gaps between top of constraint and other bounds
        # check if it is a closed contour with top and bottom files
        if (np.size(bound['file']) == 2):
            bottom_file = bound['file'][0]
            top_file = bound['file'][1]
            
            f_bottom = get_interp(bound, bottom_file, invisible=True)
            f_top = get_interp(bound, top_file, invisible=True)
            
            ax.fill_between(x, f_bottom(x), f_top(x), zorder=-1, alpha=0.3)
            line, = ax.plot(x, f_bottom(x), label=fr"{key.split('_')[0]}", alpha=1)
            ax.plot(x, f_top(x), alpha=1, color=line.get_color())


        else:
            this_file = bound['file']

            m4, Usqr = np.genfromtxt(this_file, unpack=True)

            # make sure the data points are ordered in m4
            order = np.argsort(m4)
            m4 = np.array(m4)[order]
            Usqr = np.array(Usqr)[order]

            f = get_interp(bound, this_file, invisible=True)
        
            ax.plot(x, f(x), label=fr"{key.split('_')[0]}")


    ax.set_xlim(np.min(x),np.max(x))
    ax.set_ylim(1e-8,0.5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel(fr"$|U_{4}|^2$")
    ax.set_xlabel(r"$m_{N}$/MeV")
    ax.set_yticks([ 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3,1e-2,1e-1])
    ax.grid(axis='y', which='both',dashes=(6,1),alpha=0.5,c='black',lw=0.2)

    return ax




USQR=get_combined_bounds(muon_bounds)
USQR_inv=get_combined_bounds(muon_bounds, invisible=True)

if __name__ == '__main__':

    ### Test for a few values
    ###############################
    # THIS IS THE FUNCTION THAT TAKE m4 AS INPUT
    # AND OUTPUTS THE CONSTRAINT ON UMU4^2.
    m4 = np.logspace(-3,2,1000)
    np.savetxt('umu4_constraint_new.dat', np.array([m4, USQR(m4)]).T , header='m4 (GeV) Umu4SQR ')


