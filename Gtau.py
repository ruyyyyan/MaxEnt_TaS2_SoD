import commands
import os
import sys
import time

import matplotlib
matplotlib.use('Agg')

import shutil
import os

import h5py

from pylab import *
import numpy as np
import scipy
from scipy import interpolate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

import json
import pylab
import numpy

from pylab import *
from scipy import *

from scipy          import optimize
from scipy.optimize import curve_fit

from numpy.linalg import inv

####################################################
#  parameters
####################################################
#Ms = ['bo-','rs-','g^-','mv-','c<-','y>-','kp-','bo--','rs--','g^--','mv--','c<--','y>--','kp--']

Udd = "6"
Upp = "0"
Nc = "6"
beta = "12.0"
mu = "2.2"

dens  =['5']
denss  =['5.0']

Qs = ['0']#,'0']

Ms = ['bo-','rs-','g^-','mv-','c<-','kp-','yh-']
Mss= ['bo--','rs--','g^--','mv--','c<--','y>--','kp--']

clf()
#a = loadtxt("./Gr0_U"+Udd+"_Up"+Upp+"_be"+beta+"_s1234567_N"+Nc+"_mu"+mu)
a = np.genfromtxt("./Gr0_U"+Udd+"_Up"+Upp+"_be"+beta+"_s1234567_N"+Nc+"_mu"+mu, \
                          skip_header=1, skip_footer=367)#, dtype=float)
b = np.genfromtxt("./Gr0_U"+Udd+"_Up"+Upp+"_be"+beta+"_s1234567_N"+Nc+"_mu"+mu, \
                          skip_header=123, skip_footer=245)#, dtype=float)
c = np.genfromtxt("./Gr0_U"+Udd+"_Up"+Upp+"_be"+beta+"_s1234567_N"+Nc+"_mu"+mu, \
                          skip_header=245, skip_footer=123)#, dtype=float)

plot(a[:,0], a[:,1], 'b-', linewidth = 1.5, label='Cu')
plot(b[:,0], b[:,1], 'r-', linewidth = 1.5, label='O_px')
plot(c[:,0], c[:,1], 'g-', linewidth = 1.5, label='O_py')

title("Ud="+Udd+"_Up="+Upp+"_be="+beta+"_s1234567_N="+Nc+"_mu="+mu)
grid('on')
legend(loc=1, fontsize=9.5, framealpha=1.0)
#xlim([-10,5])
#ylim([0,1.2])
xticks(fontsize=12)
yticks(fontsize=12)
xlabel(r'$\tau$',fontsize=17)
ylabel(r'$G(\tau)$',fontsize=17)
savefig("Gtau_Ud"+Udd+"Up"+Upp+"Nc"+Nc+".pdf")