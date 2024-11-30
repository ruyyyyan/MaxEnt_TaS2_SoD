import matplotlib
matplotlib.use('Agg')

import subprocess
import os
import sys
import time
import shutil

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import math

from pylab import *
import scipy
from scipy import *
from scipy import interpolate
from numpy import *

import re
import mmap
from linecache import getline

# use the common py analysis toolbox valid for all DQMC projects
cmd = "cp ../tools/dqmc_analysis_tools.py ."
os.system(cmd)

import dqmc_analysis_tools as dqmc

####################################################
#  parameters
####################################################
seed = 1234567

eh = -0.37
ets = [0.0]
tH = 0.15
tT = 0.0
Vints = [0.03, 0.1, 0.3, 0.5, 1.0]
Vints = [0.03]
Uh = 0.0
Ut = 0.1
localVs = [0.0]#, 0.5, 1.0, 1.5, 2.0, 3.0]
#localVs = [0.0]#, 3.0]


cellL = 1  # linear size in unit cell
cellN = cellL*cellL
Ncell = 4     # linear size of supercell
NL = cellL*Ncell   # linear size of whole lattice
N = NL*NL     
N3 = N*3

if cellL==1:
    Gsites = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
elif cellL==8:
    Gsites = [86, 108, 109, 111]

#betas   = [0.4, 0.45,  0.5, 0.55,  0.6,  0.7,  0.8,  0.9,   1.0, 1.2,   1.5,   1.6,   1.8,    2.0,   2.2,  2.5,  \
#           2.8,  3.0,  3.2,  3.5,  4.0,  4.2,  5.0,  6.0,   7.0, 8.0,   10.0,  12.0,  20.0]
#betas   = [4.0,  5.0,  6.0, 7.0, 8.0, 10.0, 12.0, 20.0]
betas = [6.0, 7.0, 8.0, 10.0, 12.5, 15.0]
Ls    = [60,  70,  80,  100,  100,  120]
betas = [10.0]
Ls    = [100]

Ms = ['bo-','rs-','g^-','mv-','c<-','kp-','yh-','bo--','rs--','g^--','mv--','c<--','y>--','kp--']

#########################################################################################
clf()   
T = zeros(len(betas))
batch_str = ""

cmd = 'rm Gdata/*'
os.system(cmd)

# extract various data from fortran's output files 
for Vint in Vints:
    for et in ets:
        for i in range(len(betas)):     
            #print 'beta=', betas[i]
            T[i] = 1./betas[i]
            
            for localV in localVs:
                # mus[0] for d=1 and mus[1] for 1/8 hole doping
                if cellL==1 and et==0.0:
                    if betas[i]==10.0:
                        if localV==0.0:
                            mus = [0.0, -0.2, -0.4, -0.6, -0.8, -1.0, -1.2, -1.4, -1.6]
                            mus = [0.0]
                        elif localV==0.5:
                            mus = [0.23, 1.06]
                        elif localV==1.0:
                            mus = [0.13, 0.823]
                        elif localV==1.5:
                            mus = [-0.1, 0.45]
                        elif localV==3.0:
                            mus = [-0.891, -0.75]
                    elif Vint==0.1 and betas[i]==8.0:
                        if localV==0.0:
                            mus = [0.88]
                            
                elif cellL==1 and et==1.3:
                    if Vint==9 and betas[i]==10.0:
                        if localV==0.0:
                            mus = [0.32]#, 1.24]
                        elif localV==0.5:
                            mus = [0.21, 1.13]
                        elif localV==1.0:
                            mus = [0.11, 0.88]
                        elif localV==1.5:
                            mus = [-0.14, 0.49]
                        elif localV==2.0:
                            mus = [-0.37, 0.05]
                
                # undoped and 1/8 hole doped case
                for mu in mus:
                    data = '../dqmc_SoD_test'  \
                            +'/Gr0_V'+str(Vint)+'_Uh'+str(Uh)+'_Ut'+str(Ut)+'_tH'+str(tH)+'_tT'+str(tT) \
                            +'_eh'+str(eh)+'_et'+str(et)+'_N'+str(N)+'_be'+str(betas[i]) \
                            +'_s1234567_mu'+str(mu)
                    print("Data file:", data)
                    if not os.path.isfile(data):
                        print(f"File {data} does not exist!")

                    if os.path.isfile(data):
                        #print(data)
                        Nlines = len(open(data).readlines())
                        #print Nlines 

                        # see https://www.numpy.org/devdocs/user/basics.io.genfromtxt.html for genfromtxt usage
                        #a = np.genfromtxt('../../data/'+filename, skip_header=1, skip_footer=Nlines-(Nsites2+1), \
                        #              usecols=(0,1,2,3), dtype=float)
                        #a = loadtxt('../../data/'+filename,skiprows=1)
                        #assert(len(a)==Nsites2)

                        for j in range(len(Gsites)):
                            fname = 'local_G_localV'+str(localV)+'_V'+str(Vint)+'_Uh'+str(Uh)+'_Ut'+str(Ut)+'_tH'+str(tH)+'_tT'+str(tT) \
                            +'_N'+str(N)+'_be'+str(betas[i]) \
                            +'_s1234567_mu'+str(mu)+'_site'+str(Gsites[j])

                            if os.path.isfile('./Gdata/'+fname):
                                os.remove('./Gdata/'+fname)

                            #print 'imp site ', idx_droplet[j]
                            if Gsites[j]<10:
                                phrase = str(Gsites[j])+'   '+str(Gsites[j])
                            elif Gsites[j]<100:
                                phrase = str(Gsites[j])+'  '+str(Gsites[j])
                            else:
                                phrase = str(Gsites[j])+' '+str(Gsites[j])
                            phrase = 'Gfun      '+phrase

                            # write L for maxent
                            ff = open('./Gdata/'+fname,'w') 
                            ff.write(str(Ls[i])+'\n')
                            ff.close()   # close is important to avoid abnormal issues of the file output

                            # find G(tau) lines for specific site, +2 because the Gr0 output format
                            ff = open(data, 'r')
                            dqmc.find_lines(phrase, ff, 2, Ls[i], './Gdata/'+fname)

                            if os.path.isfile('./Gdata/'+fname):
                                print (fname, 'is generated!')
