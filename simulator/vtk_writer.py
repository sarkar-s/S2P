#!/usr/bin/env python
import os
"""This file creates the vtk format lattice network plots for animation
"""
casename = []

def set_case_name(case_name):
    casename.append(case_name)

def create_vtk_template(lattice_size):
    template_file = open(casename[0]+'_template.vtk','w')
    
    print >> template_file, '# vtk DataFile Version 4.1'
    print >> template_file, 'This is a BFM network file created using SimPoP'
    print >> template_file, 'ASCII'
    print >> template_file, 'DATASET STRUCTURED_GRID'
    print >> template_file, 'DIMENSIONS '+str(lattice_size+1)+' '+str(lattice_size+1)+' '+str(lattice_size+1)
    print >> template_file, '\n'
    
    # Print the points data
    print >> template_file, 'POINTS '+str((lattice_size+1)**3)+' float'
    
    for x in range(0,lattice_size+1):
        for y in range(0,lattice_size+1):
            for z in range(0,lattice_size+1):
                print >> template_file, str(x)+' '+str(y)+' '+str(z)
    
    print >> template_file, '\n'
    
    print >> template_file, 'CELL_DATA ',((lattice_size)**3)
    print >> template_file, 'SCALARS SiteType int 1'
    print >> template_file, 'LOOKUP_TABLE default'
    
    template_file.close()
    
def create_vtk_bond_file(lattice_size,tag):
    template_file = open(casename[0]+tag+'_bonds.vtk','w')
    
    print >> template_file, '# vtk DataFile Version 1.0'
    print >> template_file, 'This is a BFM network file created using SimPoP'
    print >> template_file, 'ASCII'
    print >> template_file, 'DATASET POLYDATA'
    #print >> template_file, 'DIMENSIONS '+str(lattice_size+1)+' '+str(lattice_size+1)+' '+str(lattice_size+1)
    print >> template_file, '\n'
    
    # Print the points data
    print >> template_file, 'POINTS '+str((lattice_size+1)**3)+' float'
    
    for x in range(0,lattice_size+1):
        for y in range(0,lattice_size+1):
            for z in range(0,lattice_size+1):
                print >> template_file, str(x)+' '+str(y)+' '+str(z)
    
    template_file.close()
    
    return casename[0]+tag+'_bonds.vtk'

def case_and_step_file(tag):
    template_file = casename[0]+'_template.vtk'
    this_step_file = casename[0]+'_'+tag+'.vtk'
    
    os.system('cp '+template_file+' '+this_step_file)
    os.system('chmod 0755 '+this_step_file)
    
    return this_step_file