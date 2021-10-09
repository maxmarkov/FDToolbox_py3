import sys
sys.path.insert(0, './')

from generate_inputs import FDTcalculations
from dielectric import compute_dielectric_tensor

"""
Example showing how to compute the dielectric tensor for Eu4Sb2O
"""

##  List of computaional parameters:
# 'trans':   Threshold for translational (accoustic) modes used when inverting force constant matrix (eV/angstroem**2)
# 'rotat':   Threshold for rotational modes used when inverting elastic constant matrix (GPa)
# 'nocache': Reread the data and overwrite '_pickled.dat' cache file."
# 'asr':     Force acoustic sum rule on the FC matrix.
options = {'trans': None, 'rotat': None, 'nocache': False, 'asr': False, 'decomp': False}

# folder with calculations
folder = '/globalscratch/ucl/naps/mmarkov/Eu4Sb2O_test'

# 
fdt = FDTcalculations(filename='POSCAR', disp_ion=0.005, disp_lattice=0.00)

# compute the dielectric function
compute_dielectric_tensor(folder, fdt, options)
