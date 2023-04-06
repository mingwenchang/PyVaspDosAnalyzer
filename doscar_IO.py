#!/usr/bin/env python3
# coding=UTF-8
"""
Name: doscar_IO
Created on Fri Sep 6 10:57:09 2019
Developer: Ming-Wen Chang
E-mail: ming.wen.c@gmail.com
"""

import re
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.interpolate import make_interp_spline
from collections import OrderedDict


class Utilities:

    TDOS_channels =  OrderedDict([
                        (3, ['energy', 'tdos', 'int_tdos']),
                        (5, ['energy', 'tdos_up', 'tdos_down', 'int_tdos_up', 'int_tdos_down'])
                        ])      
    
    
    
    PDOS_channels = OrderedDict([
                    (4, ['energy','s', 'p', 'd']),  
                     
                    (5, ['energy','s', 'p', 'd', 'f']),  
                     
                    (10, ['energy','s', 'p_y', 'p_z', 'p_x', 'd_xy', 'd_yz', 'd_z^2', 'd_xz', 'd_x^2-y^2']),  
                     
                    (17, ['energy','s', 'p_y', 'p_z', 'p_x', 'd_xy', 'd_yz', 'd_z^2', 'd_xz', 'd_x^2-y^2', 
                         'f_y(3x^2-y^2)', 'f_xyz', 'f_yz^2', 'f_z^3', 'f_xz^2', 'f_z(x^2-y^2)', 'f_x(x^2-3y^2)']),
                    
                    (7, ['energy','s_up', 's_down', 'p_up', 'p_down', 'd_up', 'd_down']),
                    
                    (9, ['energy','s_up', 's_down', 'p_up', 'p_down', 'd_up', 'd_down', 'f_up', 'f_down']),
                    
                    (19, ['energy','s_up', 's_down', 'p_y_up', 'p_y_down', 'p_z_up', 'p_z_down', 'p_x_up', 'p_x_down', 
                        'd_xy_up', 'd_xy_down', 'd_yz_up', 'd_yz_down', 'd_z^2_up', 'd_z^2_down', 'd_xz_up', 
                        'd_xz_down', 'd_x^2-y^2_up', 'd_x^2-y^2_down']), 
                    
                    (33, ['energy','s_up', 's_down', 'p_y_up', 'p_y_down', 'p_z_up', 'p_z_down', 'p_x_up', 'p_x_down',
                        'd_xy_up', 'd_xy_down', 'd_yz_up', 'd_yz_down', 'd_z^2_up', 'd_z^2_down', 'd_xz_up',
                        'd_xz_down', 'd_x^2-y^2_up', 'd_x^2-y^2_down', 
                        'f_y(3x^2-y^2)_up', 'f_y(3x^2-y^2)_down', 'f_xyz_up', 'f_xyz_down', 'f_yz^2_up', 'f_yz^2_down', 
                        'f_z^3_up', 'f_z^3_down', 'f_xz^2_up', 'f_xz^2_down', 'f_z(x^2-y^2)_up', 'f_z(x^2-y^2)_down', 
                        'f_x(x^2-3y^2)_up', 'f_x(x^2-3y^2)_down'])
                    ])   
 
    # up to 24 colors stored   
    try:
        import seaborn as sns
        sns.set_style("white")
        sns.set_context("paper") #,font_scale=0.5, rc={"lines.linewidth": 20, "grid.linewidth": 0.00})
        clist = sns.color_palette("muted", 24)
    except:
        clist =  ['blue', 'red', 'gold', 'salmon', 'lightcoral', 'lightskyblue', 'darkgreen', 'black',
                 'orange','powderblue','olivedrab', 'burlywood',  'indianred', 
                'steelblue', 'lawngreen', 'y', 'hotpink', 'slategrey', 'yellowgreen','palegreen', 
                'sandybrown', 'tomato', 'darkviolet', 'lightgreen', 'tan','maroon']

    cdict = OrderedDict([
            ('s', 'slategrey'), ('p', 'salmon'), ('d', 'gold'), ('f','lightgreen'), ('tdos', 'powderblue')
            ])

                    
    def set_spin_polarized_channel(self, names):  
        #names = list(names)
        sp_channel = [ ]
        for name in names:
            sp_channel.append(name +'_up')
            sp_channel.append(name +'_down')
            
        if  'energy' in sp_channel[0]:
            sp_channel.remove('energy_up')
            sp_channel.remove('energy_down')
            sp_channel.insert(0, 'energy')
        return sp_channel
    
    def reduce_to_undecomposed_channel(self, names):
        channel = [ ]
        for name in names:
            if name == 'energy':
                channel.append(name)
            else:
                name = name.split('_')[0] + '_' + name.split('_')[-1]
                if name not in channel :
                    channel.append(name)
        return channel


class Doscar:
    
    def __init__(self, filename='DOSCAR', ispin=None, lmax=None, lorbit=None, 
                 eshift=0.00, undecomposed=True, mode='vasp'):  
        
        self.filename = filename
        self.ispin = ispin 
        self.lmax = lmax
        self.lorbit = lorbit
        self.eshift = eshift  
        self.undecomposed = undecomposed
        self.analdbc = False
        self.tdos = None
        self.pdos = None
        
        if os.path.exists(self.filename):
            if mode.lower()[0] == 'v':
                self.read_vasp_DOSCAR()
            elif mode.lower()[0]  == 'l':
                self.read_lobster_DOSCAR()
            else:
                raise ValueError("vasp and lobster modes are supported")
                
        else:
            raise IOError("%s doesn't exist" %(self.filename))
                

    def read_vasp_DOSCAR(self, filename=None):
        uti = Utilities()
        if filename is None:
            filename = self.filename 
        nheadlines = 6    
        self.data = []
        with open(filename, 'r' ) as txt: 
            #Read basic information
            lines = [txt.readline() for i in range(nheadlines)]
            self.natoms = int(lines[0].split()[0]) #numbar of atoms in this system
            self.nedos = int(lines[5].split()[2])  #number of gridpoints on which the DOS is evaluated
            self.efermi = float(lines[5].split()[3]) #fermi level
            self.emax = float(lines[5].split()[0]) - self.efermi
            self.emin = float(lines[5].split()[1]) - self.efermi
            
            #Read tdos information
            tdos = []
            for j in range(self.nedos):
                values = [float(value) for value in txt.readline().split()]
                tdos.append(values)
                
            channel = len(tdos[0])
            names = uti.TDOS_channels[channel] 
            df = pd.DataFrame(data=tdos, columns=names)
            
            

            self.energies = df['energy'] - self.efermi
            self.efermi -= self.efermi
            
            if self.eshift > 0.00 : 
                self.shift_fermi_level(self.eshift)
            
            df = self.invert_dos_values_of_spin_down(df)
            self.data.append(df)
            self.tdos = self.data[0]

            #Read atomic pdos  
            for i in range(1, self.natoms+1):
                txt.readline().split()
                pdos = [ ]
                for j in range(self.nedos):
                    values = [float(value) for value in txt.readline().split()]
                    pdos.append(values)
                    
                channel = len(pdos[0])
                names = uti.PDOS_channels[channel] 
                df = pd.DataFrame(data=pdos, columns=names)
                
                if self.undecomposed:
                    df = self.reduce_to_undecomposed(df)

                df['energy'] = self.energies
                df = self.invert_dos_values_of_spin_down(df)                
                self.data.append(df)
 
    def read_lobster_DOSCAR(self, filename=None):
        uti = Utilities()
        if filename is None:
            filename = self.filename 
        nheadlines = 6    
        self.data = []
        with open(filename, 'r' ) as txt: 
            #Read basic information
            lines = [txt.readline() for i in range(nheadlines)]
            self.natoms = int(lines[0].split()[0]) #numbar of atoms in this system
            self.nedos = int(lines[5].split()[2])  #number of gridpoints on which the DOS is evaluated
            self.efermi = float(lines[5].split()[3]) #fermi level
            self.emax = float(lines[5].split()[0]) - self.efermi
            self.emin = float(lines[5].split()[1]) - self.efermi
            
            #Read tdos information
            tdos = []
            for j in range(self.nedos):
                values = [float(value) for value in txt.readline().split()]
                tdos.append(values)
                
            channel = len(tdos[0])
            names = uti.TDOS_channels[channel] 
            df = pd.DataFrame(data=tdos, columns=names)
            
            self.energies = df['energy'] #- self.efermi
            self.efermi -= self.efermi
            
            if self.eshift > 0.00: 
                self.shift_fermi_level(self.eshift)
                
            df = self.invert_dos_values_of_spin_down(df)
            self.data.append(df)
            self.tdos = self.data[0]
            
            #Read atomic pdos  
            for i in range(1, self.natoms+1):
                info = txt.readline().split()
                names = info[7:]
                names.insert(0, 'energy')

                pdos = [ ]
                for j in range(self.nedos):
                    values = [float(value) for value in txt.readline().split()]
                    pdos.append(values)

                if len(pdos[0]) - 1 >= len(names): 
                    names = uti.set_spin_polarized_channel(names)
                
                df = pd.DataFrame(data=pdos, columns=names)
                if self.undecomposed:
                    df = self.reduce_to_undecomposed(df)
                    
                df = self.invert_dos_values_of_spin_down(df)
                df['energy'] = self.energies
                self.data.append(df)
                
               
    def reduce_to_undecomposed(self, df0):
        uti = Utilities()
        names = uti.reduce_to_undecomposed_channel(df0.columns)
        df = pd.DataFrame(data=0.00, index=range(self.nedos), columns=names)
        
        for name in names[1:]:
            for column in df0.columns:
                channel = column.split('_')[0] + '_' + column.split('_')[-1]
                if channel == name:
                    df[name] += df0[column]
        df['energy'] = self.energies         
        return df
    
    def shift_fermi_level(self, eshift):
        self.efermi = eshift 
        self.energies -= self.efermi  
        self.emin = self.energies[0]
        self.emax = self.energies[-1]
        print ('Fermi level is shifted to %s eV' %(eshift))
        print ('All other energies are also shifted with respect to the new Fermi level')
    
    def invert_dos_values_of_spin_down(self, df):
        for column in df.columns:
            if 'down' in column:
                df[column] *= -1
        return df

    def get_tdos(self):
        self.tdos = self.data[0]
        return self.tdos
    
    def get_atomic_pdos(self, idx=0):
        return self.data[idx] 

    def get_pdos(self, kwargs):        
        frames = [ ]
        for key in kwargs.keys():
            atomtype = key
            orbital = kwargs[key][0] 
            atomic_list = kwargs[key][1]
            df = self.sum_atomic_pdos(atomic_list=atomic_list, data=self.data)
            if orbital is not None:
                df = self.select_orbital(df=df, name=orbital)
            del df['energy']
            #Change names of columns
            df.columns = [atomtype + '-' + name for name in df.columns]
            frames.append(df)
        pdos = pd.concat(frames, axis=1)   
        pdos.insert(0, 'energy', self.energies)
        self.pdos = pdos
        return self.pdos
   
    def sum_atomic_pdos(self, atomic_list=[0], data=None):
        if data is None:
            data = self.data
        #initialize a zero datafram
        num = atomic_list[0]
        names = data[num].columns
        df = pd.DataFrame(data=0.00, index=range(self.nedos), columns=names)
        for num in atomic_list:
            df += data[num].iloc[:,1:] 
        df['energy'] =  self.energies
        return df
        
    def select_orbital(self, df=None, name='d', mode='fuzzy'):
        if df is None:
            df = self.pdos
            
        ids = [0]
        
        for idx, key in enumerate(df.columns):
            if mode[0].lower() == 'e': #explicit
                if key == name:
                    ids.append(idx)
            elif mode[0].lower() == 'f': #fuzzy
                key = key.split('_')[0].split('-')
                if len(key) > 1:
                    key = key[-1]
                else:
                    key = key[0]
                    
                key = re.sub('\d', '', key)
                if key == name:
                    ids.append(idx)
            else: #normal 
                key = key.split('_')[0] 
                if key == name:
                    ids.append(idx)

        df = df.iloc[:,ids]
        return df        
    
    def split_doscar(self):
        for idx, df in enumerate(self.data):
            filename='dos' + str(idx)
            self.save_df_to_txt(df, filename)
        
    def save_dos(self, filename='pdos'):
        if self.pdos is None:
            df = self.get_tdos()
            filename = 'tdos'
        else:
            df = self.pdos
        self.save_df_to_txt(df, filename)
    
    def save_df_to_txt(self, df, filename='df.txt'):
        with open(filename, 'w') as txt:
            txt.write(df.to_string(index=False))
             
    def calculate_dbc(self, df=None, orbital='d', erange=None):
        if df is None:
            df = self.pdos
        
        if erange is None:
            erange = (self.emin, self.efermi+2.00)
        mask = (self.energies >= erange[0]) & (self.energies <= erange[1])
            
        df = self.select_orbital(df, orbital)
        
        x = self.energies[mask]
        if len(df.keys()) == 3:
            y1 = df.iloc[:,1].values
            y2 = df.iloc[:,2].values
        else:# len(df.keys()) == 2:
            y1 = df.iloc[:,1].values
            y2 = df.iloc[:,1].values
  
        y1 = y1[mask]
        y2 = y2[mask]
        
        nele_up = simps(y1, x) 
        nele_down = simps(y2, x)
        self.nelectrons_in_pdos =  abs(nele_up) , abs(nele_down)
        
        dbc_up = simps(y1*x, x) / nele_up 
        dbc_down = simps(y2*x, x) / nele_down 
        self.dbc =  dbc_up, dbc_down
        self.analdbc = True    
        return self.dbc    

   
    def plot_pdos(self, df=None, filename='pdos.png', erange=None, smooth=True,
                  size=(1.5, 6), dpi=150, cmap=None, line_width=0.05, line_alpha=0.50, 
                  fill=True, fill_alpha=0.65, fontsize=10, style='column'):

        if df is None:
            df = self.pdos
        
        if erange is None:
            emin, emax  = (self.efermi-6.00, self.efermi+4.00)
        else:
            emin, emax = erange[0], erange[1] 
            
        if cmap is None:
            cmap = Utilities().clist
        coloridx = 0

        plt.figure(figsize=size, dpi=dpi)
        for key in df.columns[1:]:
            if 'int' in key : 
                continue

            #Color setting
            label = key.split('_')[0]
            if type(cmap) is dict: 
                try:
                    cc = cmap[label]
                except:
                    orb = re.sub("\d", "", label.split('-')[-1])
                    cc = cmap[orb] 
            else:
                if 'down' not in key:
                    cc = cmap[coloridx]
                else:
                    coloridx -= 1
                    cc = cmap[coloridx]
                coloridx += 1
 
            #Values for plotting
            evalues = self.energies  
            dosvalues = df[key].values
            if smooth:
                evalues, dosvalues = self.smooth_line(evalues, dosvalues, grid=self.nedos*10)
            
            mask = (evalues >= emin) & (evalues <= emax)            
            evalues= evalues[mask]
            dosvalues = dosvalues[mask] 
            
            dmax = 1.05 * np.max(np.abs(dosvalues))  
            dmin = -dmax
            
            #style = 'd'
            if style.lower()[0] == 'c':
                x = dosvalues
                y = evalues
            else:
                x = evalues
                y = dosvalues
                
            if 'down' not in key:  
                plt.plot(x, y, linewidth=line_width, alpha=line_alpha, label=label, color=cc)
                legend =plt.legend(bbox_to_anchor=(1.00, 1.02), loc= 'upper left', prop={'size':fontsize*0.8}, frameon=True)
                [i.set_linewidth(3) for i in legend.legendHandles]
            else:
                plt.plot(x, y, linewidth=line_width, alpha=line_alpha, label='', color=cc)
                
            if fill: 
                if style.lower()[0] == 'c':
                    plt.fill_betweenx(y, x, -0.000,  where=y <= 0, interpolate=True, color=cc, alpha=fill_alpha) 
                    plt.fill_betweenx(y, x, +0.000,  where=y >= 0, interpolate=True, color=cc, alpha=fill_alpha*0.3)
                else:
                    plt.fill_between(x, y, -0.000,  where=x <= 0, interpolate=True, color=cc, alpha=fill_alpha) 
                    plt.fill_between(x, y, -0.000,  where=x >= 0, interpolate=True, color=cc, alpha=fill_alpha*0.3) 
                
        
        if style.lower()[0] == 'c':
            #Plot data in the column style
            plt.ylim([emin, emax])
            plt.xlim([dmin, dmax])
            plt.ylabel( '$E - E_{f}\ (eV)$', size=fontsize)
            plt.xlabel( 'pDOS (a.u.)', size=fontsize)
            plt.tick_params(which='major', labelleft= True, left=True, labelbottom=False, direction ='out', labelsize=fontsize) 
            plt.axvline(x=0, linewidth=1.00, linestyle='-', color='gray', alpha=0.80)
            plt.axhline(y=self.efermi, linewidth=2.50, linestyle='-', color='pink', alpha=0.80)
            if self.analdbc:
                plt.axhline(y=max(self.dbc), xmin=0.00, xmax=1.00, linewidth=2.00, linestyle='-', color='green', alpha=0.50)
        else:
            #Plot data in the row style
            plt.xlim([emin, emax])
            plt.ylim([dmin, dmax])
            plt.xlabel( '$E - E_{f}\ (eV)$', size=fontsize)
            plt.ylabel( 'pDOS (a.u.)', size=fontsize)
            plt.tick_params(which='major', labelleft=False,  direction ='in', labelsize=fontsize) 
            plt.axvline(x=self.efermi, linewidth=2.50, linestyle='-', color='pink', alpha=0.80)
            plt.axhline(y=0, linewidth=0.20, linestyle='-', color='gray', alpha=0.90)
            if self.analdbc:
                plt.axvline(x=max(self.dbc), ymin=0.00, ymax=1.00, linewidth=2.00, linestyle='-', color='green', alpha=0.50)

        plt.savefig(filename, bbox_inches="tight")
                         
        
    def smooth_line(self, x, y, grid=None):
        if grid is None:
            grid = len(x) * 100
        xnew = np.linspace(x.min(), x.max(), grid) 
        bsplobj = make_interp_spline(x, y, k=3) #BSpline object
        ynew = bsplobj(xnew) #smoothed data
        return xnew, ynew

debug = True 
if debug:
    atoms = list(range(1, 8+1))
    kwargs = {'Pt':('d', [5])}
    a = Doscar('DOSCAR', mode='v')
    a.get_pdos(kwargs)
    a.calculate_dbc(erange=(-6,2))
    a.plot_pdos(erange=(-6,2))
    a.save_dos()
    print(a.dbc)
    print(a.nelectrons_in_pdos)

