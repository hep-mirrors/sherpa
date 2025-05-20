#Input the random seed ending * that identifies the run, from runcard_*.yaml
#Must be string and 3 characters, include 00 for 001, 002 etc.
random_seed_end = "001"


#Import functions and set design parameters
import numpy as np 
from matplotlib import pyplot as plt
from matplotlib import rc
import glob as glob

from functions import *

plt.rcParams["axes.grid"] = False

rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

#### Design properties ####
fScaling = 1.5 #Decides overall scaling of figures, 1 should be optimal
fWidth = fScaling*4.7747
fHeight = fScaling*3.81976
fontSize = 10.5*fScaling
fontSizeLegend = 9.5*fScaling

path_end = "_s*" + random_seed_end + ".dat"

path_start_2D = "../HIAnalysis_2D/"

histo_name = "t+mg2"
x, y, y_err, y_err_2, N_groups = AverageFrom3DHisto(path_start_2D + histo_name + path_end)

if N_groups > 0:
    plt.figure(figsize=(fWidth,fHeight))
    plt.title(r'$m_g^2$ vs t', fontsize=fontSize)


    x = x*GeV_to_fm
    y = y/GeV_to_fm
    y_err = y_err*GeV_to_fm
    plt.plot(x, y, color='tab:blue', linewidth=1.5, label=r'$\hat{q}$')
    errorband_bottom = y-y_err
    errorband_bottom[errorband_bottom<0] = 0.
    plt.fill_between(x, errorband_bottom, y+y_err, alpha=0.5, edgecolor='tab:blue', facecolor='tab:blue', label=r'Error mean between ' + str(N_groups) + ' event groups')

    plt.ylabel(r'$m_g^2$ [$\mathrm{GeV}^2$]')
    plt.xlabel(r'$t$ [1/\mathrm{GeV}]', fontsize=fontSize)
    plt.legend(loc="upper right",prop={"size":fontSizeLegend})
    plt.tight_layout()
else:
    print("Could not find any data for histograms: " + histo_name)