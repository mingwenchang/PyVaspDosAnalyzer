#!/usr/bin/env python3
# coding=UTF-8
"""
Name: pypltpdos
Created on Fri Sep 6 10:57:09 2019
Developer: Ming-Wen Chang
E-mail: ming.wen.c@gmail.com
"""

import re
import sys
#import vasp_modules.doscar_IO as dio
from collections import OrderedDict
import testdos as dio


class State:
    #DOSCAR parameters
    dosfile = 'DOSCAR' 
    system = ''
    efermi = None
    undecomposed = True
    splitdos = False
    savepdos = True
    kwargs = OrderedDict() # {}
    
    #DBC 
    dbc = False
    dbcorbital= 'd'
    dbcintegrange = None
    
    #Plot parameters 
    erange = None
    smooth = True
    plotpdos = True
    figuresize = (6, 1.5)
    dpi = 150
    cmap = None
    width1 = 0.50
    alpha1 = 1.00
    fill = True
    alpha2 = 0.50
    fontsize = 12
    neat = False
    stack='row'
    
            
def analyze_atom_object(string):
    import re
    string = re.sub('type:|orbital:|atomlist:|\s', '', string)
    split_string = re.split(',',string)
    
    atomtype = split_string[0]
    
    atomorbital = split_string[1].lower()
    if atomorbital == 'all' or atomorbital == 'none':
        atomorbital = None
    
    atomlist = []
    for term in split_string[2:]:
        if '-' in term:
            start = int(term.split('-')[0])
            end = int(term.split('-')[-1])
            [atomlist.append(_) for _ in range(start, end+1)]
        else:
            term = int(term)
            atomlist.append(term)
    #remove duplicates in  atomlist       
    atomlist = list(dict.fromkeys(atomlist))
    atomtype, atomorbital, atomlist
    return atomtype, atomorbital, atomlist
    
     
def read_input(filename='ini.dos'):
    with open(filename, 'r') as txt:
        for line in txt:
            #Skip lines start with: '#', '!', '\s', '\n'
            skip = any([line.startswith(_) for _ in ['#', '!', '\s', '\n']])
            if not skip:
                line = re.sub('(?<!\'|\")#.*', '', line)
                line = re.sub('\s|\'|\"', '', line)
                line = re.split('=', line)
                keyword = line[0].upper()
                arguments = line[1]
                print(keyword, ':', arguments)
                
                if keyword == 'SYSTEM':
                    if arguments != '':
                        state.jobname = arguments + '-' 
                        
                elif keyword == 'FILENAME':
                    state.dosfile = arguments
                                        
                elif keyword == 'EFERMI':
                    state.efermi = float(arguments)
                    
                elif keyword == 'ATOM_OBJECT':
                    atomtype, atomorbital, atomlist = analyze_atom_object(arguments)
                    state.kwargs[atomtype] = (atomorbital, atomlist)
                    
                elif keyword == 'UNDECOMPOSED': 
                    if arguments.capitalize() == 'True':
                        state.undecomposed = True
                    else:
                        state.undecomposed = False  
                        
                elif  keyword == 'SPLITDOS':
                    if arguments.capitalize() == 'True':
                        state.splitdos = True  
                    else:
                        state.splitdos = False 
                        
                elif  keyword == 'SAVEPDOS':
                    if arguments.capitalize() == 'True':
                        state.savepdos = True
                    else:
                        state.savepdos = False

                elif keyword == 'ANALDBC':
                    if arguments.capitalize() == 'True':
                        state.dbc = True 
                    else:
                        state.dbc = False
                        
                elif keyword == 'DBCORBITAL':
                    state.dbcorbital = arguments.lower()

                elif keyword == 'DBCINTEGRANGE':
                    arguments = arguments.strip(')').strip('(').split(',')
                    if arguments[0].lower()[0] == 'd':
                        state.dbcintegrange  = None
                    else:
                        state.dbcintegrange = (float(arguments[0]), float(arguments[1]))
                        
                elif  keyword == 'PLOTPDOS':
                    if arguments.capitalize() == 'True':
                         state.plotpdos = True
                    else:
                        state.plotpdos = False
                        
                elif keyword == 'SMOOTH':
                    if arguments.capitalize() == 'True':
                         state.smooth = True
                    else:
                        state.smooth = False
                        
                elif keyword == 'FILLAREA': 
                    if arguments.capitalize() == 'True':
                        state.fill = True
                    else:
                        state.fill = False
                        
                elif  keyword == 'NEAT':
                    if arguments.capitalize() == 'True':
                         state.neat = True
                    else:
                        state.neat = False
                        
                elif keyword == 'STACK':
                    state.stack = arguments.lower()
                    
                elif keyword == 'PLOTERANGE':
                    arguments = arguments.strip(')').strip('(').split(',')
                    if arguments[0].lower()[0] == 'd':
                        state.erange =  None
                    else:
                        state.erange = (float(arguments[0]), float(arguments[1]))
                        
                elif  keyword == 'FIGURESIZE':
                    arguments = arguments.strip(')').strip('(').split(',')
                    state.figuresize = (float(arguments[0]), float(arguments[1]))
                    
                elif keyword == 'LINEWIDTH':
                    state.width1 = float(arguments)
                
                elif keyword == 'FILLALPHA':
                    state.alpha2 = float(arguments)
                    
                elif keyword == 'DPI':
                    state.dpi = float(arguments)
                    
                elif keyword == 'FONTSIZE':
                    state.fontsize = float(arguments)
                    
                elif keyword == 'CMAP':
                    if arguments[0].lower()[0] == 'd':
                        state.cmap =  None
                    else:
                        cmap = [term.strip(',').strip('[').strip(']') for term in arguments.split(',')]
                        state.cmap = cmap    
                else:
                    raise ValueError('KEYWORD: %s is not supported' %(keyword))
            
# Main Function   
if __name__ == "__main__":
    
    if len(sys.argv) > 1:
        inp = sys.argv[1]
    else:
        inp = 'pdos.ini'
    
    
    print("o=======================================================o")
    print("            Reading DOSCAR and pdos.ini files            ")
    print("o=======================================================o")
    
    state = State()
    read_input(filename=inp)

    doscar = dio.Doscar(filename=state.dosfile, efermi=state.efermi, 
                        undecomposed=state.undecomposed,)
    
    if state.splitdos:
        print("o=======================================================o")
        print("                 Splitting pDOS                          ")
        print("o=======================================================o")
        
        doscar.split_doscar()
        print("DOSCAR file has been splitted")    
    
    doscar.get_pdos(state.kwargs)
    if state.dbc:
        print("o=======================================================o")
        print("                 Analyzing %s-band center                " %(state.dbcorbital))
        print("o=======================================================o")
        doscar.calculate_dbc(orbital=state.dbcorbital, erange=state.dbcintegrange) 
        dbc = doscar.dbc
        nelectrons = sum(doscar.nelectrons_in_pdos)
        with open(state.jobname+'dbc.dat', 'w') as txt:
            txt.write("%s-band center: %s eV,  %s eV\n" %(state.dbcorbital, dbc[0], dbc[1]))
            txt.write("There are %s electrons filling the %s-band" %(nelectrons, state.dbcorbital))
        print ("%s-band center: %s eV,  %s eV" %(state.dbcorbital, dbc[0], dbc[1]))
        print("There are %s electrons filling the %s-band" %(nelectrons, state.dbcorbital))
        
    if state.savepdos:
        doscar.save_dos(doscar.pdos, state.jobname+'pdos.dat')
        print("Your PDOS data have been saved as %s file" %(state.jobname+'pdos.dat'))
        
    if state.plotpdos:
        print("o=======================================================o")
        print("                 Plotting pDOS                          ")
        print("o=======================================================o")
        doscar.plot_pdos(filename=state.jobname+'pdos.png', erange=state.erange, 
                         smooth=state.smooth, size=state.figuresize, 
                         dpi=state.dpi, cmap=state.cmap, 
                         line_width=state.width1, line_alpha=state.alpha1 ,
                         fill=state.fill, fill_alpha=state.alpha2, fontsize=state.fontsize,
                         stack=state.stack, neat=state.neat)
        print('pDOS plot has been saved as %s' %(state.jobname+'pdos.png'))


        

         

    
    
