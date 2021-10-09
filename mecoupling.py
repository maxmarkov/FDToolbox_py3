#!/usr/bin/env python

from fdtoolbox.calculation_set import calculation_set
from generate_inputs import FDTcalculations
from fdtoolbox.utility import set_common_options, mat2str
from fdtoolbox.linear_expansion import linear_expansion

import sys
import os
import optparse
import numpy as np

def compute_magnetoelectric_tensor(folder, fdt, saxis, options):
    '''
    Compute the magnetoelectric tensor from the set of magnetization & polarization calculations.
    folder: folder with VASP calculations.
    fdt: object
    saxis: quantization axis as a list
    options: dictionary with computational parameters
    '''

    usecache = False
     
    print (f'Using saxis: {saxis}')
    print (f'Using cache: {usecache}')
    print(f'Reading from directory: {folder}')
    
    cs=calculation_set.read_directory(folder, fdt, saxis, usecache, 'nscf')
    
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
    
    me, units = lin_exp.magneto_electric_coupling(lattice=True)
    print('Full ME coupling {units}')
    print(mat2str(me))
