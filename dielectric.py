#!/usr/bin/env python

from fdtoolbox.calculation_set import calculation_set
from generate_inputs import FDTcalculations 
from fdtoolbox.utility import set_common_options, add_common_options, mat2str
from fdtoolbox.linear_expansion import linear_expansion
import sys
import os
import optparse

saxis = [0, 0, 1]
usecache = False

options = {'trans': None, 'rotat': None, 'nocache': False, 'asr': False, 'decomp': False}
arguments = []

directory = '/globalscratch/ucl/naps/mmarkov/Eu4Sb2O_test'

fdt = FDTcalculations(filename='/CECI/home/ucl/modl/mmarkov/tutorial_FDToolbox/reference_initial_files/POSCAR', disp_ion=0.005, disp_lattice=0.00)

print ('Using saxis: %s' % saxis)
print ('Using cache: %s' % usecache)

print('Reading from directory: %s' % directory)

cs=calculation_set.read_directory(directory, fdt, saxis, usecache, 'scf')

cs.set_ionic('.*/calc_.._.*[0123]')
cs.set_lattice('.*/calc_00_...')
cs.set_groundstate('.*/calc_00_0')

calculation_set.BPH_MAXDIST = 0.25
calculation_set.BPH_CORRECTION = 0.5

cs.try_fix_displacements()  
cs.try_fix_polarizations()
lin_exp = linear_expansion(cs)

# ADD OPTIONS
set_common_options(cs, lin_exp, options)

lin_exp.calculate_expansion_coefficients(update_from_calcset=True)

chitens, units = lin_exp.electric_susceptibility(lattice = False)
print(f'Dielectric tensor ionic contribution {units}')
print(mat2str(chitens))

chitens, units = lin_exp.electric_susceptibility()
print(f'Dielectric tensor ionic and cell contribution {units}')
print(mat2str(chitens))
