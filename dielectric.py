#!/usr/bin/env python

from fdtoolbox.calculation_set import calculation_set
from generate_inputs import FDTcalculations 
from fdtoolbox.utility import set_common_options, mat2str
from fdtoolbox.linear_expansion import linear_expansion
import sys
import os
import optparse


def compute_dielectric_tensor(folder, fdt, options):
    '''
    Compute the dielectric tensor from the set of polarization calculations.
    folder: folder with VASP calculations.
    fdt: object
    options: dictionary with computational parameters
    '''

    usecache = False

    print (f'Using cache: {usecache}')
    print(f'Reading from directory: {folder}')

    cs=calculation_set.read_directory(folder, fdt, [0, 0, 1], usecache, 'scf')

    cs.set_ionic('.*/calc_.._.*[0123]')
    cs.set_lattice('.*/calc_00_...')
    cs.set_groundstate('.*/calc_00_0')

    calculation_set.BPH_MAXDIST = 0.25
    calculation_set.BPH_CORRECTION = 0.5

    cs.try_fix_displacements()  
    cs.try_fix_polarizations()
    lin_exp = linear_expansion(cs)

    # ADD OPTIONS
    if len(options) == 0:
       # 'trans':   Threshold for translational (accoustic) modes used when inverting force constant matrix (eV/angstroem**2)
       # 'rotat':   Threshold for rotational modes used when inverting elastic constant matrix (GPa)
       # 'nocache': Reread the data and overwrite '_pickled.dat' cache file."
       # 'asr':     Force acoustic sum rule on the FC matrix.
       options = {'trans': None, 'rotat': None, 'nocache': False, 'asr': False, 'decomp': False}
    set_common_options(cs, lin_exp, options)

    lin_exp.calculate_expansion_coefficients(update_from_calcset=True)

    chitens, units = lin_exp.electric_susceptibility(lattice=False)
    print(f'\nDielectric tensor ionic contribution {units}')
    print(mat2str(chitens))

    chitens, units = lin_exp.electric_susceptibility(lattice=True)
    print(f'\nDielectric tensor ionic and cell contribution {units}')
    print(mat2str(chitens))
