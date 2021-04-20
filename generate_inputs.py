#!/usr/bin/env python

import copy
import sys
import os
import argparse
import shutil

from pymatgen.io.vasp.sets import MPStaticSet
from pymatgen.transformations.site_transformations import TranslateSitesTransformation
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from tools.submit import generate_submit_manneback


def move_all_atoms(structure, disp_ion):
    """
    Create structures with displaced ions
    """
    yield 'calc_00_0', structure
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
    """
    Create structures with displaced lattice
    """
    for axnum, axname in enumerate(['a', 'b', 'c']):
      for dirnum, dirname in enumerate(['X', 'Y', 'Z']):
        for orientation in ['+', '-']:
           structure_displaced = copy.deepcopy(structure)
           structure_displaced.lattice.matrix[axnum, dirnum] += disp_lattice
           yield 'calc_00_%s%s%s' % (axname, orientation, dirname), structure_displaced
           disp_lattice = -disp_lattice

def generate_scf_folders(structures, user_incar):
    """
    Generate calculations from the structures dictionary 
    """
    for name in structures:
        if os.path.exists(name)==False:
           os.mkdir(name)
        structures[name].to(fmt='POSCAR', filename='/'.join([name, 'POSCAR']))
        MPStaticSet(structures[name], user_incar_settings=user_incar, user_kpoints_settings={"reciprocal_density": 500}).write_input(name)
        print(name+'\n')

def generate_old_polarization_folders(structures):
    """
    Generate the polarization calculations (old format) from the completed scf calculation
    """
    for name in structures:
        for i in range(3):
            # compute the Berry phase along the x-i, y- and z-directions (IGPAR = 1,2,3)
            # Note: calculations usualy are very sensetive to the choice of "NPPSTR" parameter. You might consider to reduce it.
            MPStaticSet.from_prev_calc(name, user_incar_settings={"LBERRY":True,"IGPAR":i+1,"ICHARG":2,"NPAR":2,"ALGO":"Fast","NPPSTR":6,"LAECHG":False,"LVHAR":False, "ISIF": 2}).write_input(os.path.join(name,'Berry_'+str(i+1)))
            #shutil.copy(os.path.join(name,'CHGCAR'), os.path.join(name,'Berry_'+str(i+1)))
            shutil.copy(os.path.join(name,'KPOINTS'), os.path.join(name,'Berry_'+str(i+1)))


class FDTcalculations(object):
    """
    This class generates the self-consistent calculations for FDToolbox 
    """
    def __init__(self, filename='POSCAR', disp_ion=0.005, disp_lattice=0.000):
       """
       Initializes the class from the original structure file and displacement amplitudes
       """
       self.filename = filename
       self.disp_ion = disp_ion
       self.disp_lattice = disp_lattice 
       if os.path.isfile(filename) == False:
          sys.exit('Correct path to a structure POSCAR/cif file must be provided')
    
       self.structure = Structure.from_file(self.filename)
       print(self.structure)

       # dictionary with ionic displacements
       self.structures_ion = {}
       if self.disp_ion > 0.0:
          for name, structure_displaced in move_all_atoms(self.structure, self.disp_ion):
              self.structures_ion[name] = structure_displaced
       else:
          print("Ionic displacement is zero. Only the reference structure is returned.\n")

       # dictionary with lattice displacements
       self.structures_lattice = {}
       if self.disp_lattice > 0.0:
          for name, structure_displaced in move_all_lattice(self.structure, self.disp_lattice):
              self.structures_lattice[name] = structure_displaced
       else:
          print("Lattice displacement is zero. Strain contribution won't be computed.\n")

    def generate_scf(self, user_incar={"ALGO":"Normal","EDIFF":5e-8,"LWAVE":False,"LVHAR":False,"LAECHG":False,"NPAR":4,"ICHARG":2}):
        """
        Generate folders with calculations: INCAR, POSCAR, POTCAR, KPOINTS
        """
        if len(self.structures_ion) > 0:
           generate_scf_folders(self.structures_ion, user_incar)
        if len(self.structures_lattice) > 0:
           generate_scf_folders(self.structures_lattice, user_incar)

    def check_scf(self):
        """
        Check either all scf calculations have been completed correctly 
        """
        calcs_failed = []
        if len(self.structures_ion) > 0:
           for folder in self.structures_ion:
               f = os.path.join(folder,'CHGCAR')
               if os.path.exists(f) == False or os.path.getsize(f) < 1000:
                  calcs_failed.append(folder)
                  print('Calculation in {} folder has not been completed'.format(folder))
        if len(self.structures_lattice) > 0:
           for folder in self.structures_lattice:
               f = os.path.join(folder,'CHGCAR')
               if os.path.exists(f) == False or os.path.getsize(f) < 1000:
                  calcs_failed.append(folder)
        if len(calcs_failed) == 0:
           return(True)
        else:
           print('Following calculations have problems: ', calcs_failed)
           return(False)


    def generate_old_polarization(self):
        """
        Generate folders with (old-type) polarization calculations for each direction
        """
        if len(self.structures_ion) > 0:
           generate_old_polarization_folders(self.structures_ion)
        if len(self.structures_lattice) > 0:
           generate_old_polarization_folders(self.structures_lattice)


####         STEP 1           ###
# customize the default settings for INCAR file
#user_incar_custom = {"ALGO":"Normal","LDAUU":{"Eu":6.0},"LDAUL":{"Eu":3, "Sb":-1, "O":-1},"EDIFF":5e-8, "LWAVE":False,"LVHAR":False,"LAECHG":False,"LREAL":False,"NPAR":4,"ICHARG":2}

#fdt = FDTcalculations(filename='/home/ucl/modl/mmarkov/soft/FDToolbox_py3_dev/TEST_FOLDER/POSCAR', disp_ion=0.005, disp_lattice=0.00)
#fdt.generate_scf(user_incar=user_incar_custom)

#folders = list(fdt.structures_ion.keys())
#generate_submit_manneback(folders, nprocs = 16, path_to_exec="/home/ucl/modl/mmarkov/soft/vasp/vasp.5.4.1_couplageEM/bin/vasp_std")


####         STEP 2           ###
fdt = FDTcalculations(filename='/home/ucl/modl/mmarkov/soft/FDToolbox_py3_dev/TEST_FOLDER/POSCAR', disp_ion=0.005, disp_lattice=0.00)

scf_completed = fdt.check_scf()
print('Scf calculations completed succesfully? ', scf_completed)

#if scf_completed:
#   fdt.generate_polarization_old()
