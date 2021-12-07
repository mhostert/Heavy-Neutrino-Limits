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
# Which constraints do you want to exclude from the final combination
EXCLUDE_THESE_CONSTRAINTS = []


PATH = 'digitized'
files_umu4 = [
        [PATH+"/umu4/BEBC.dat", PATH+"/Umu4/BEBC_top.dat"],
        [PATH+"/umu4/CHARMII.dat", PATH+"/Umu4/CHARMII_top.dat"],
        PATH+"/umu4/DELPHI_1997.dat",
        [PATH+"/umu4/FMMF.dat", PATH+"/Umu4/FMMF_top.dat"],
        PATH+"/umu4/K_to_numu.dat",
        PATH+"/umu4/L3.dat",
        [PATH+"/umu4/NA3.dat",PATH+"/Umu4/NA3_top.dat"],
        [PATH+"/umu4/NuTeV.dat", PATH+"/Umu4/NuTeV_top.dat"],
        [PATH+"/umu4/PS191_bottom.dat", PATH+"/Umu4/PS191_top.dat"],
        PATH+"/peak_searches/E949.dat",
        PATH+"/peak_searches/KEK_K_nu_mu.dat",
        PATH+"/peak_searches/NA62_(2017)_K_nu_mu.dat",
        PATH+"/peak_searches/NA62_2018_K_nu_mu.dat",
        PATH+"/deGouvea/lowmassUmu4.dat",
        PATH+"/Kusenko/Umu4.dat",
        PATH+"/Kusenko/Umu4_v2.dat",
        PATH+"/Kusenko/Umu4_2.dat",
        PATH+"/Kusenko/Umu4_3.dat",
        PATH+"/Kusenko/Umu4_4.dat",
        PATH+"/Kusenko/Umu4_5.dat",
        PATH+"/umu4/SuperK_2019.dat",
        PATH+"/umu4/T2K_2019.dat",
        PATH+"/umu4/NA62_2021.dat",
        PATH+"/umu4/PIENU_2019_lowT.dat",
        PATH+"/umu4/PIENU_2019_highT.dat",
        PATH+"/umu4/PSI_1987.dat",
        PATH+"/umu4/DELPHI_short_95CL.dat",
        PATH+"/umu4/DELPHI_long_95CL.dat",
        PATH+"/umu4/SIN_1987.dat",
        ]


## Return a function that takes m4 [GeV] and gives bound on Ualpha4^2
def func_bound(listoffiles):
    GREYCOLOR = "grey"
    ALP = 0.5
    y=[]

    x=np.logspace(-3,2,1000)
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



        if i in [9,10,11,12,14,15,16,17,18,19,21,23,24,25,28]:
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




## Return a function that takes m4 [GeV] and gives bound on Ualpha4^2
def individual_bounds(ax, listoffiles):

    x=np.logspace(-3,2,1000)
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



        if i in [9,10,11,12,14,15,16,17,18,19,21,23,24,25,28]:
            ax.plot(m4*1e-3, Umu4sq, c=colors[i], ls =linestyles[i], lw=2.0, label=labels[i])
            f = interpolate.interp1d(m4*1e-3, Umu4sq, kind='linear', bounds_error=False, fill_value=1.0, assume_sorted=False)    
        else:
            ax.plot(m4, Umu4sq, c=colors[i], ls =linestyles[i], lw=2.0, label=labels[i])
            f = interpolate.interp1d(m4, Umu4sq, kind='linear', bounds_error=False, fill_value=1.0, assume_sorted=False)    


    return ax


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
            "red"]*3

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
            r"deGouvea"]*3

linestyles = ["-","-","-","-","-","-","-","-","-",":",":",":",":",":"]*3

USQR=func_bound(files_umu4)

if __name__ == '__main__':

    ### Test for a few values
    ###############################
    # THIS IS THE FUNCTION THAT TAKE m4 AS INPUT
    # AND OUTPUTS THE CONSTRAINT ON UMU4^2.
    m4 = np.logspace(-3,2,1000)
    np.savetxt('umu4_constraint_v1.dat', np.array([m4, USQR(m4)]).T , header='m4 (GeV) Umu4SQR ')




#### Some old code for sensitivity lines

# def plot_sensitivities(axis, LS=(5,0.01), SHIP=1, FCC=1, MATHUSLA=1, SBN=1,NA62=1, DUNE=1, SBNCOLOR='purple', DUNECOLOR='blue'):

#     LW = 1.6
#     # SHiP
#     if SHIP ==1:
#         m4, Umu4sq = np.genfromtxt(PATH+'/1805.08567/lower_SHIP.dat', unpack=True)
#         m4_up, Umu4sq_up = np.genfromtxt(PATH+'/1805.08567/upper_SHIP.dat', unpack=True)
#         axis.plot(m4,Umu4sq, color='black', lw = LW, zorder=300, dashes=LS)
#         axis.plot(m4_up,Umu4sq_up, color='black', lw = LW, zorder=300, dashes=LS)

#     if MATHUSLA ==1:
#         # MATHUSLA
#         m4, Umu4sq = np.genfromtxt(PATH+'/1805.08567/mathusla_lower.dat', unpack=True)
#         m4_up, Umu4sq_up = np.genfromtxt(PATH+'/1805.08567/mathusla_upper.dat', unpack=True)
#         axis.plot(m4,Umu4sq, color='darkorange', lw = LW, zorder=300, dashes=LS)
#         axis.plot(m4_up,Umu4sq_up, color='darkorange', lw = LW, zorder=300, dashes=LS)

#     if FCC ==1:
#         # FCC-ee
#         m4, Umu4sq = np.genfromtxt(PATH+'/1805.08567/FCCee_lower.dat', unpack=True)
#         m4_up, Umu4sq_up = np.genfromtxt(PATH+'/1805.08567/FCCee_upper.dat', unpack=True)
#         axis.plot(m4,Umu4sq, color='darkgreen', lw = LW, zorder=300, dashes=LS)
#         axis.plot(m4_up,Umu4sq_up, color='darkgreen', lw = LW, zorder=300, dashes=LS)

#     if SBN ==1:
#         m4, Umu4sq = np.genfromtxt(PATH+'/tommaso_nuphys/Umu4_SBN.dat', unpack=True)
#         axis.plot(m4,Umu4sq, color=SBNCOLOR, lw = 2.0, zorder=300, dashes=LS)

#     if NA62 ==1:
#         m4, Umu4sq = np.genfromtxt(PATH+'/NA62/NA62.dat', unpack=True)
#         m4_up, Umu4sq_up = np.genfromtxt(PATH+'/NA62/NA62_upper.dat', unpack=True)
#         axis.plot(m4,Umu4sq, color='red', lw = LW, zorder=300, dashes=LS)
#         axis.plot(m4_up,Umu4sq_up, color='red', lw = LW, zorder=300, dashes=LS)

#     if DUNE ==1:
# #         m4, Umu4sq = np.genfromtxt('../digitized/tommaso_nuphys/DUNE_both.dat', unpack=True)
# #         axis.plot(m4,Umu4sq, color='darkblue', lw = LW, zorder=300)
#         m4, Umu4sq = np.genfromtxt(PATH+'/tommaso_DUNE/Umu4.dat', unpack=True)
#         m4_up, Umu4sq_up = np.genfromtxt(PATH+'/tommaso_DUNE/Umu4_upper.dat', unpack=True)
#         axis.plot(m4,Umu4sq, color=DUNECOLOR, lw = LW, zorder=300, dashes=LS)
#         axis.plot(m4_up,Umu4sq_up, color=DUNECOLOR, lw = LW, zorder=300, dashes=LS)

