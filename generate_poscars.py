#!/usr/bin/env python

import numpy as np

import copy
import sys
import os
import shutil

sys.path.append(os.path.join(sys.path[0],'../'))
from fdtoolbox.calculation_set import calculation
from pymatgen.io.vasp.sets import MPStaticSet

disp_length = [0.005, 0.005, 0.005]
lat_disp_length = [0.01, 0.01, 0.01]

def move_all_atoms(calc):
  yield 'calc_00_0', calc
  for numat, at in enumerate(calc.atoms):
    for direction in range(3):
      displacement = np.zeros((1,3))
      displacement[0, direction] = disp_length[direction % len(disp_length)]
      
      for orientation in ['+', '-']:
        c = copy.deepcopy(calc)
        c.atoms[numat] += np.squeeze(displacement)
        yield 'calc_%02d_%s%d' % (numat+1, orientation, direction+1), c
      
        displacement *= -1

def move_all_lattice(calc):
  for axnum, axname in enumerate(['a', 'b', 'c']):
    for dirnum, dirname in enumerate(['X', 'Y', 'Z']):
      displacement = np.zeros((1,3))
      displacement[0, dirnum] = lat_disp_length[dirnum % len(lat_disp_length)]
      for orientation in ['+', '-']:
        c = copy.deepcopy(calc)
        c.unit_cell[axnum] += np.squeeze(displacement)
        c.atoms = np.dot( calc.atoms, np.dot( np.linalg.inv( calc.unit_cell ), c.unit_cell ) )
        yield 'calc_00_%s%s%s' % (axname, orientation, dirname), c
      
        displacement *= -1

def generate_all_calcs(calc):
  for name, c in move_all_atoms(calc):
    yield name, c
  for name, c in move_all_lattice(calc):
    yield name, c

# sys.argv is a list, which contains the command-line arguments passed to the script.
print(sys.argv)
if len(sys.argv) < 2:
  name = (sys.argv[0].split('/'))[-1]
  msg = ('{sep}Usage:{sep}'
         '{x} filename ion_displacement_lengths lattice_displacement_lengths{sep}'
         '  filename - name of the POSCAR/CONTCAR file containing optimised system{sep}'
         '  ion_displacement_lengths - ion displacement in angstroms single value for uniform displacement, or quoted triplet of values{sep}'
         '  lattice_displacement_lengths - same as ion_displacement_lengths but for lattice vectors{sep}{sep}'
         'Example:{sep}'
         '{x} ./POSCAR 0.005 "0.01 0.01 0.05"{sep}'.format(x=name, sep='\n'))
  print(msg) 
  sys.exit()

system = calculation()
system.load_from_poscar(sys.argv[1])


if len(sys.argv) > 2:
  disp_length = [ float(d) for d in sys.argv[2].split() ]

if len(sys.argv) > 3:
  lat_disp_length = [ float(d) for d in sys.argv[3].split() ]

generator = generate_all_calcs
if max(disp_length) == 0.:
  generator = move_all_lattice
if max(lat_disp_length) == 0.:
  generator = move_all_atoms

for name, calc in generator(system):
  try:
    os.mkdir(name)
  except:
    pass
  calc.save_to_poscar('/'.join( [name, 'POSCAR'] ))
#  MPStaticSet(calc.structure, user_incar_settings={"ALGO":"Normal","LDAUU":{"Eu":6.0},"LDAUL":{"Eu":3},"EDIFF":5e-8,"LWAVE":False,"LVHAR":False,"LAECHG":False,"NPAR":4,"ICHARG":2}).write_input(name)
  print(name)
