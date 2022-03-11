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

##
# Which constraints to exclude from the final combination
EXCLUDE_THESE_CONSTRAINTS = []


PATH = os.path.dirname(os.path.abspath(__file__))+'/digitized'
files_utau4 = [ 
        PATH+"/utau4/Borexino_RPlestid.dat",
        PATH+"/utau4/NOMAD_UL_90CL_2002.dat",
        PATH+"/utau4/CHARM_UL_90CL_2002.dat",
        PATH+"/utau4/T2K_UL_90CL_2019.dat",
        PATH+"/utau4/DELPHI_short_95CL.dat",
        PATH+"/utau4/DELPHI_long_95CL.dat",
        PATH+"/utau4/EWPD_Bolton.dat",
        PATH+"/utau4/Bdecays_Bolton.dat",
        PATH+"/utau4/l_universality_Bolton.dat",
        PATH+"/utau4/DELPHI_1997.dat",
        ]
        
## Return a function that takes m4 [GeV] and gives bound on Umu4^2

def func_bound(listoffiles):
    GREYCOLOR = "grey"
    ALP = 0.5
    y=[]

    x=np.logspace(-3,3,1000)
    for i in range(np.size(listoffiles)):


        ## 
        # EXCLUDE CERTAIN CONSTRAINTS THAT DO NOT APPLY
        filename = listoffiles[i]
        if (np.size(filename) == 2):
            name = filename[0]
        else:
            name = filename

        m4, Umu4sq = np.genfromtxt(name, unpack=True)

        order = np.argsort(m4)
        m4 = np.array(m4)[order]
        Umu4sq = np.array(Umu4sq)[order]



        if i in [0,1,2,3]:
#       axis.plot(m4*1e-3, Umu4sq, c=colors[i], ls =linestyles[i], lw=2.0, label=labels[i])
            f = interpolate.interp1d(m4*1e-3, Umu4sq, kind='linear', bounds_error=False, fill_value=1.0, assume_sorted=False)    
            y.append(f(x))
        else:
#             axis.plot(m4, Umu4sq, c=colors[i], ls =linestyles[i], lw=2.0, label=labels[i])
            f = interpolate.interp1d(m4, Umu4sq, kind='linear', bounds_error=False, fill_value=1.0, assume_sorted=False)    
            y.append(f(x))

    y = np.array(y)
    z = np.ones(np.size(y[0]))
    for i in range(0,np.size(listoffiles)):
        for j in range(0, np.size(y[i])):
            if i not in EXCLUDE_THESE_CONSTRAINTS:
                if y[i,j] < z[j]:
                    z[j] = y[i,j]

    return interpolate.interp1d(x, z, kind='linear', bounds_error=False, fill_value=1.0, assume_sorted=False)    

    
###############################
# THIS PLOTS ALL CONSTRAINTS -- JUST AN APPLICATION OF THE FUNCTION ABOVE.
colors = [  "black",
            "indigo",
            "lightblue",
            "orange",
            "red",
            "lightgreen",
            "pink",
            "darkblue",
            "darkgreen",
            "black",
            "indigo",
            "orange",
            "red",
            "red"]

labels = [  r"BEBC",
            r"CHARMII",
            r"DELPHI",
            r"FMMF",
            r"$K \to \nu \mu$",
            r"L3",
            r"NA3",
            r"NuTeV",
            r"PS191",
            r"E949",
            r"KEK",
            r"NA62 (2017)",
            r"NA62 (2018)",
            r"deGouvea"]

linestyles = ["-","-","-","-","-","-","-","-","-",":",":",":",":",":"]

USQR=func_bound(files_utau4)
if __name__ == '__main__':
        


    ### Test for a few values
    ###############################
    # THIS IS THE FUNCTION THAT TAKE m4 AS INPUT
    # AND OUTPUTS THE CONSTRAINT ON UMU4^2.
    m4 = np.logspace(-3,3,1000)
    np.savetxt('utau4_constraint_v1.dat', np.array([m4, USQR(m4)]).T , header='m4 (GeV) Utau4SQR ')