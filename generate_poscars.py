#!/usr/bin/env python

import copy
import sys
import os
import numpy as np
import argparse

#from pymatgen.io.vasp.sets import MPStaticSet
from pymatgen.transformations.site_transformations import TranslateSitesTransformation
from pymatgen.core.structure import Structure

def move_all_atoms(structure, disp_ion):
  yield 'calc_00_0', structure
  print(structure.lattice)
  if disp_ion == 0.:
     print("Ionic displacement is zero. Only the reference structure is returned.\n") 
     return
  atoms = [i for i in range(len(structure))]
  for atom in atoms:
    for direction in range(3):
      displacement = [0.]*3
      displacement[direction] = disp_ion
      
      for orientation in ['+', '-']:
         transform = TranslateSitesTransformation([atom], displacement, vector_in_frac_coords=True)
         structure_displaced = transform.apply_transformation(structure)
         yield 'calc_%02d_%s%d' % (atom+1, orientation, direction+1), structure_displaced      
         displacement[direction] = -displacement[direction]

def move_all_lattice(structure, disp_lattice):
  if disp_lattice == 0.:
     print("Lattice displacement is zero. Strain contribution won't be computed.\n")
     return
  for axnum, axname in enumerate(['a', 'b', 'c']):
    for dirnum, dirname in enumerate(['X', 'Y', 'Z']):
      lattice_matrix = copy.deepcopy(structure.lattice.matrix)
      lattice_matrix[axnum, dirnum] += disp_lattice
      print(axname, dirname)
      print(lattice_matrix)
#      for orientation in ['+', '-']:
#        c = copy.deepcopy(calc)
#        c.unit_cell[axnum] += np.squeeze(displacement)
#        c.atoms = np.dot( calc.atoms, np.dot( np.linalg.inv( calc.unit_cell ), c.unit_cell ) )
#        yield 'calc_00_%s%s%s' % (axname, orientation, dirname), c
#      
#        displacement[dirnum] = -displacement[dirnum]

def generate_all_calcs(structure, disp_ion, disp_lattice):
  for name, c in move_all_atoms(structure, disp_ion):
    yield name, c
  move_all_lattice(structure, disp_lattice)
#  for name, c in move_all_lattice(structure, disp_lattice):
#    yield name, c

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--filename', type=str, default='POSCAR', help='File with original structure. POSCAR or cif formats')
    parser.add_argument('--disp_ion', type=float, default='0.005', help='Magnitude of ionic displacement')
    parser.add_argument('--disp_lattice', type=float, default='0.000', help='Magnitude of lattice vectors displacements')
    opt = parser.parse_args()

    # customize the default settings for INCAR file
    user_incar_custom = {"ALGO":"Normal","LDAUU":{"Eu":6.0},"LDAUL":{"Eu":3},"EDIFF":5e-8,"LWAVE":False,"LVHAR":False,"LAECHG":False,"NPAR":4,"ICHARG":2} 

    if os.path.isfile(opt.filename) == False:
       sys.exit('Correct path to a structure POSCAR/cif file must be provided after --filename')
    
    structure = Structure.from_file(opt.filename)

    for name, structure_displaced in generate_all_calcs(structure, opt.disp_ion, opt.disp_lattice):
        try:
          os.mkdir(name)
        except:
          pass
        structure.to(fmt='POSCAR', filename='/'.join( [name, 'POSCAR']))
        print(name)
        #MPStaticSet(structure, user_incar_custom).write_input(name)
