import sys
sys.path.insert(0, './')

import os
from generate_inputs import FDTcalculations
from tools.submit import generate_submit_manneback

"""
Computational STEP 2:
        Example showing how to generate folders for the polarization (noSOC) and magnetization (SOC) calculations
"""

#
fdt = FDTcalculations(filename='/home/ucl/modl/mmarkov/soft/FDToolbox_py3_dev/TEST_FOLDER/POSCAR', disp_ion=0.005, disp_lattice=0.00)

scf_completed = fdt.check_scf()
print('All SCF calculations completed succesfully? ', scf_completed)

# customize the default settings for INCAR file: provide the Hubbard model parameters for each element
user_incar_nscf_custom = {"LDAUU":{"Eu":6, "Sb": 0, "O": 0},"LDAUL":{"Eu":3, "Sb":-1, "O":-1}}

# quabtization axis
saxis = (0,0,1)

# path to the root folder of VASP
path_bin = "/home/ucl/modl/mmarkov/soft/vasp/vasp.5.4.1_couplageEM/"

if scf_completed:
   fdt.generate_old_polarization()
   fdt.generate_new_polarization()
   fdt.generate_nscf_soc(saxis, user_incar_nscf_custom)

   folders = list(fdt.structures_ion.keys()) + list(fdt.structures_lattice.keys())
   for subfolder in ['Berry_1','Berry_2','Berry_3','Berry_new']:
       generate_submit_manneback(folders, subfolder=subfolder, nprocs = 16, path_to_exec=os.path.join(path_bin,"bin/vasp_std"))
   generate_submit_manneback(folders, subfolder='nscf_SOC', nprocs = 16, path_to_exec=os.path.join(path_bin,"bin/vasp_ncl"))
