# -*- coding: utf-8 -*-
"""
Created on Fri May 10 09:26:06 2019

@author: 20180239
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


names = ['dz_10.00', 'dz_3.50', 'dz_3.00', 'dz_2.50', 'dz_2.12']
lables=[ 'energy', 'up', 'down' ]

#Set Figure size
fig, axs = plt.subplots(nrows=1, ncols=len(names),  sharey=True, figsize=(len(names)*3, 6), dpi=200)
fig.subplots_adjust(wspace=0.05)
#Set  x-y range
xrange = [-3.00, 3.00]
yrange = [-8.00, 3.00]


#read data
for i, name in enumerate(names):
    filename = 'C:\\Users\\20180239\\Desktop\\COHP-plot\\COHP\\plane\\%s' %(name)
    
    #Read data 
    df = pd.read_csv( filename, skiprows=5, header = None, usecols=[0,1,3], delim_whitespace=True)
    df.columns = lables
    
    ax =axs[i]
    x1   = -1*df.up 
    x2   = -1*df.down
    y1   = df.energy
    y2   = df.energy
     

    if i == 0:
        text_dz = r'$10.0\ \AA $'
        x = x1
        y = y1
        xy= (0.67, 0.95)
        cc1 = 'red' #'red' or 'cyan'
        cc2 = 'pink' # 'pink' or 'lightblue'
        xycoords = 'axes fraction'
    elif i == 1:
        text_dz = r'$3.50\ \AA $'
        x = x1
        y = y1
    elif i == 2:
        text_dz = r'$3.00\ \AA $'
        x = x1
        y = y1
    elif i == 3:
        text_dz = r'$2.50\ \AA $'
        x = x1
        y = y1
    elif i == 4:
        text_dz = r'$2.12\ \AA $'
        x = x1
        y = y1
    elif i == 5:
        text_dz = r'$2.12\ \AA $'
        x = x1
        y = y1
        
    """
    else:  
        text_dz ="  GM  "
        xy= (0.65, 0.93)
        xycoords = 'axes fraction'
    """    
    

    #Plot
    emask = (y >= yrange[0]) & (y <= yrange[1])
    xmax = np.max([abs(np.min(x)), abs(np.max(x))])
    xrange = [-xmax*1.10, xmax*1.10]
    
    ax.plot(x, y, alpha=1.00, label=None, color='gray', linewidth=0.50)
    ax.fill_betweenx(y, x, -0.005,  where=y <= 0, interpolate=True, color=cc1, alpha=0.80) 
    ax.fill_betweenx(y, x, +0.000,  where=y >= 0, interpolate=True, color=cc2, alpha=0.80) 
    
    
    """ 
    x = x2
    y = y2
    ax.plot(x, y, alpha=1.00, label=None, color='gray', linewidth=0.50)
    ax.fill_betweenx(y, x, -0.005,  where=y <= 0, interpolate=True, color='blue', alpha=0.80) 
    ax.fill_betweenx(y, x, +0.000,  where=y >= 0, interpolate=True, color='lightblue', alpha=0.80) 
    """ 

    
    #Set annotate 
    text_Ef = r'$\epsilon_F$' 
    bbox=dict(boxstyle='square', pad=0.5, fc='w', ec='none', alpha=1.00)    
    ax.annotate(text_Ef, xy=(0.750, 0.710),  xycoords=xycoords, color= 'k',
                bbox=bbox, fontsize=20)
    

    
    bbox=dict(boxstyle='square', pad=0.25, fc=cc1, ec='none', alpha=0.70)
    ax.annotate(text_dz, xy=xy, xycoords=xycoords, color= 'w',
                bbox=bbox, fontsize=15)

    #Set ticks, legend box and X-Y label
    ax.set_xlim(xrange[0], xrange[1])
    ax.set_ylim(yrange[0], yrange[1])
    ax.set_xlabel('-pCOHP', fontsize=15)
    

    ax.set_xticks(ticks=[])
    ax.annotate('anti-bonding', xy=(0.02, 0.015), xycoords='axes fraction',
                color= 'k',fontsize=10)
    
    ax.annotate('bonding', xy=(0.73, 0.015), xycoords='axes fraction',
                color= 'k',fontsize=10)
    


    #ax.set_xticklabels(['anti', 'bond'])
    #Set lines for Fermi level and y=0
    ax.axvline(x=0.00, ymin=0.00, ymax=1.00, linewidth=1.20, linestyle='-', color='black')
    ax.axhline(y=0.00, xmin=0.00, xmax=1.00, linewidth=1.20, linestyle='-', color='black') 



axs[0].set_ylabel('Energy (eV)', fontsize=15)
