#!/usr/bin/env python

from fdtoolbox.calculation_set import calculation_set
from generate_inputs import FDTcalculations 
from fdtoolbox.utility import set_common_options
from fdtoolbox.linear_expansion import linear_expansion
import sys
import os
#import optparse

saxis = [0, 0, 1]
usecache = False

directory = '/globalscratch/ucl/naps/mmarkov/Eu4Sb2O_test' #'/globalscratch/ucl/naps/mmarkov/fdtoolbox'

fdt = FDTcalculations(filename='/home/ucl/modl/mmarkov/soft/FDToolbox_py3_dev/TEST_FOLDER/POSCAR', disp_ion=0.005, disp_lattice=0.00)

print ('Using saxis: %s' % saxis)
print ('Using cache: %s' % usecache)

print('Reading from directory: %s' % directory)

cs=calculation_set.read_directory(directory, fdt, saxis, usecache)

cs.set_ionic('.*/calc_.._.*[0123]')
cs.set_lattice('.*/calc_00_...')
cs.set_groundstate('.*/calc_00_0')

calculation_set.BPH_MAXDIST = 0.25
calculation_set.BPH_CORRECTION = 0.5

cs.try_fix_displacements()  
cs.try_fix_polarizations()
lin_exp = linear_expansion( cs )

# ADD OPTIONS
#set_common_options(cs, lin_exp, options)

lin_exp.calculate_expansion_coefficients()


chitens, units = lin_exp.electric_susceptibility(lattice = False)
print('Dielectric tensor ionic contribution (%s)'%units)
print(mat2str( chitens ))

#chitens, units = lin_exp.electric_susceptibility()
#print('Dielectric tensor ionic and cell contribution (%s)'%units)
#print(mat2str( chitens ))

#if options.decomp:
##  evals = linalg.eig(lin_exp.B_m_n)[0]
#  print('enum\teval\tcontribution')
#  for i in range(lin_exp.B_m_n.shape[0]):
#    inv_B_m_n, _, ee = safe_inv( lin_exp.B_m_n, cs.TRANSLATIONAL_MODE_THRESHOLD, iterate_components = [i] )
#    print('%d\t%f\t%s'%(i, ee[i]*cs.groundstate.volume, mat2str( ( EEAEV_TO_EPSILON0*dot(lin_exp.B_m_alpha.T, dot(inv_B_m_n, lin_exp.B_m_alpha) )).flatten() )[:-1] ))
