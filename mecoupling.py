#!/usr/bin/env python

from fdtoolbox.calculation_set import calculation_set
from generate_inputs import FDTcalculations
from fdtoolbox.utility import set_common_options, add_common_options, mat2str
from fdtoolbox.linear_expansion import linear_expansion

import sys
import os
import optparse
import numpy as np

saxis = [0, 0, 1]
usecache = False
 
print ('Using saxis: %s' % saxis)
print ('Using cache: %s' % usecache)

directory = '/globalscratch/ucl/naps/mmarkov/Eu4Sb2O_test'
print('Reading from directory: %s' % directory)

fdt = FDTcalculations(filename='/CECI/home/ucl/modl/mmarkov/tutorial_FDToolbox/reference_initial_files/POSCAR', disp_ion=0.005, disp_lattice=0.00)

cs=calculation_set.read_directory(directory, fdt, saxis, usecache, 'nscf')

cs.set_ionic('.*/calc_.._.*[0123]')
cs.set_lattice('.*/calc_00_...')
cs.set_groundstate('.*/calc_00_0')

calculation_set.BPH_MAXDIST = 0.25
calculation_set.BPH_CORRECTION = 0.5

cs.try_fix_displacements()
cs.try_fix_polarizations()

lin_exp = linear_expansion(cs)
lin_exp.calculate_expansion_coefficients()

me, units = lin_exp.magneto_electric_coupling(lattice=False)
print(f'Clamped cell ME coupling {units}')
print(mat2str(me))
#print(f'Maximum coupling: {np.sqrt(max(np.linalg.eig(np.dot(me,me.T))[0])} {units}')

me, units = lin_exp.magneto_electric_coupling()
print('Full ME coupling {units}')
print(mat2str(me))
