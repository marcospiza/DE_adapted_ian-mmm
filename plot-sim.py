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

    
x_fig = 40.0 / 2.54
y_fig = 20.0 / 2.54
font = {'size' : 16}
matplotlib.rc('font', **font)


fig = plt.figure(figsize=(x_fig,y_fig))
gs = gridspec.GridSpec(1,3,wspace=0.5,hspace=.2,width_ratios=(2,2,1) )

ax = fig.add_subplot(gs[0]) 
ax1 = fig.add_subplot(gs[1])


ax.scatter(RFA,PBobs,facecolor='none',edgecolor='C0',s=80)
ax.scatter(RFA,PBopt,marker="v",facecolor='none',edgecolor='C1',s=80)

ax.set_xlabel("RFA, MJ m$^{-2}$ d$^{-1}$")
ax.set_ylabel("PB, g CO$_2$ m$^{-2}$")

ax1.scatter(PBobs,PBopt,facecolor='none',edgecolor='C0',s=80)
ax1.set_ylabel(r"PB$_{\mathrm{opt}}$, g CO$_2$ m$^{-2}$")
ax1.set_xlabel(r"PB$_{\mathrm{obs}}$, g CO$_2$ m$^{-2}$")

# -------------Calculating correlation coefficient
corr_matrix = np.corrcoef(PBobs, PBopt)
corr = corr_matrix[0,1]
R_sq = corr**2

str_R =  "$r^2$ = " + " {:15.3f}".format(R_sq)
print(str_R)
#ax1.text()





erro = 0
for ibs,ipt in zip(PBobs,PBopt):
    erro = erro + (ibs - ipt)**2
#    print(erro)
print("------------------------------")
print( np.sqrt(erro)/len(PBobs),len(PBobs))
plt.savefig("error_comp.tiff",dpi=800,bbox_inches='tight')
