#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 17:13:51 2022

@author: marcos
"""

import matplotlib
#from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mtick



# Reading file of results and parameters

f = open("de-result.out","r")
lines = f.readlines()
f.close()

flag = True
i = 1
while(flag):
    flag = not("*****" in lines[i])
    print(lines[i])
    i = i + 1

RFA,PBobs,PBopt = [],[],[]
for l in lines[i+1:]:
    l = l.split()
    RFA.append(float(l[0]))
    PBobs.append(float(l[1]))
    PBopt.append(float(l[2]))

    
x_fig = 20.0 / 2.54
y_fig = 10.0 / 2.54
font = {'size' : 14}
matplotlib.rc('font', **font)


fig = plt.figure(figsize=(x_fig,y_fig))
gs = gridspec.GridSpec(1,3,wspace=0.5,hspace=.2,width_ratios=(2,2,1) )

ax = fig.add_subplot(gs[0]) 
ax1 = fig.add_subplot(gs[1])


ax.scatter(RFA,PBobs)
ax.scatter(RFA,PBopt)

ax1.scatter(PBopt,PBobs)

erro = 0
for ibs,ipt in zip(PBobs,PBopt):
    erro = erro + (ibs - ipt)**2
#    print(erro)
print("------------------------------")
print( np.sqrt(erro)/len(PBobs),len(PBobs))
