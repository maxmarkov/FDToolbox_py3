import sys
sys.path.insert(0, './')

from generate_inputs import FDTcalculations
from tools.submit import generate_submit_manneback

"""
Computational STEP 1:
	Example showing how to generate folders for the ground state calculations
"""

# customize the default settings for INCAR file: provide the Hubbard model parameters for each element
user_incar_scf_custom = {"LDAUU":{"Eu":6.0},"LDAUL":{"Eu":3, "Sb":-1, "O":-1}}

#
fdt = FDTcalculations(filename='POSCAR', disp_ion=0.005, disp_lattice=0.00)

# generate folders with INCAR, POSCAR, KPOINTS and POTCAR files for each displaced structure 
fdt.generate_scf(user_incar_scf_custom)

folders = list(fdt.structures_ion.keys())

# Optional: generate submit scripts for the Manneback cluster. 
generate_submit_manneback(folders, subfolder='', nprocs = 16, path_to_exec="/home/ucl/modl/mmarkov/soft/vasp/vasp.5.4.1_couplageEM/bin/vasp_std")
