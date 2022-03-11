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
EXCLUDE_THESE_CONSTRAINTS = [22]


PATH = os.path.dirname(os.path.abspath(__file__))+'/digitized'
files_ue4 = [ 
        PATH+"/ue4/NA62_2020.dat",
        PATH+"/ue4/PIENU.dat",
        PATH+"/ue4/KEK.dat",
        PATH+"/ue4/PIENU_Bolton.dat",
        PATH+"/ue4/kev_bounds/Bound1.txt",
        PATH+"/ue4/kev_bounds/Bound2.txt",
        PATH+"/ue4/kev_bounds/Bound3.txt",
        PATH+"/ue4/kev_bounds/Bound4.txt",
        PATH+"/ue4/kev_bounds/Bound5.txt",
        PATH+"/ue4/kev_bounds/Bound6.txt",
        PATH+"/ue4/kev_bounds/Bound7.txt",
        PATH+"/ue4/kev_bounds/Bound8.txt",
        PATH+"/ue4/kev_bounds/Bound9.txt",
        PATH+"/ue4/kev_bounds/Bound10.txt",
        PATH+"/ue4/kev_bounds/Bound11.txt",
        PATH+"/ue4/kev_bounds/Bound12.txt",
        PATH+"/ue4/F.dat",
        PATH+"/ue4/Fermi2_F.dat",
        PATH+"/ue4/EWPD_Bolton.dat",
        PATH+"/ue4/CHARM_1986.dat",
        PATH+"/ue4/DELPHI_1997.dat",
        PATH+"/ue4/T2K_2019_profiled.dat",
        PATH+"/ue4/TRIUMF_1992.dat",
        PATH+"/ue4/SuperK_2019.dat",
        PATH+"/ue4/KENU_bryman_shrock.dat",
        PATH+"/ue4/PIENU_bryman_shrock.dat",
        PATH+"/ue4/PIENU-H_bryman_shrock.dat",
        PATH+"/ue4/De2_bryman_shrock.dat",
        PATH+"/ue4/DELPHI_short_95CL.dat",
        PATH+"/ue4/DELPHI_long_95CL.dat",
        # PATH+"/ue4/DELPHI.dat",
        ]


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



        if i in [0,1,2,16,17,21,22,24,25,26,27]:
            # axis.plot(m4*1e-3, Umu4sq, c=colors[i], ls =linestyles[i], lw=2.0, label=labels[i])
            f = interpolate.interp1d(m4*1e-3, Umu4sq, kind='linear', bounds_error=False, fill_value=1.0, assume_sorted=False)    
            y.append(f(x))
        elif i in [4,5,6,7,8,9,10,11,12,13,14]:
            # axis.plot(m4*1e-3, Umu4sq, c=colors[i], ls =linestyles[i], lw=2.0, label=labels[i])
            f = interpolate.interp1d(m4*1e-9, Umu4sq*1e-2, kind='linear', bounds_error=False, fill_value=1.0, assume_sorted=False)    
            y.append(f(x))
        elif i in [15]:
            # axis.plot(m4*1e-3, Umu4sq, c=colors[i], ls =linestyles[i], lw=2.0, label=labels[i])
            f = interpolate.interp1d(m4*1e-6, Umu4sq*1e-2, kind='linear', bounds_error=False, fill_value=1.0, assume_sorted=False)    
            y.append(f(x))
        else:
            # axis.plot(m4, Umu4sq, c=colors[i], ls =linestyles[i], lw=2.0, label=labels[i])
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

USQR=func_bound(files_ue4)
if __name__ == '__main__':
        
    ### Test for a few values
    ###############################
    # THIS IS THE FUNCTION THAT TAKE m4 AS INPUT
    # AND OUTPUTS THE CONSTRAINT ON UMU4^2.
    m4 = np.logspace(-3,3,1000)
    np.savetxt('ue4_constraint_v1.dat', np.array([m4, USQR(m4)]).T , header='m4 (GeV) Ue4SQR ')