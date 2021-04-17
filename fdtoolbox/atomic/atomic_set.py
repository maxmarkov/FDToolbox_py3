import numpy as np
from fdtoolbox.atomic.utility import loggable
from pymatgen.core.structure import Structure

class atomic_set(loggable):
  """
  A class for keeping the information about set of atoms in a unit cell.
  - self.unit_cell - unit cell (3x3 A)
  - self.recip_cell - reciprocal unit cell (3x3 A**-1)
  - self.volume - volume calculated from the unit cell (1 A**3)
  - self.name - name of the system as put in POSCAR
  - self.num_per_type - list of species numbers 
  - self.num_atoms - number of atoms
  - self.atoms - atomic positions internal coordinates (self.num_atomsx3 1)
  """  
  #def get__volume(self):
  #  if not hasattr(self,'__volume'):
  #    self.__volume = np.linalg.det(self.unit_cell)
  #  return self.__volume    
  #volume = property( get__volume, None )

  def load_from_poscar(self, filename):
    """
    Loads the system from POSCAR using pymatgen.
    """
    structure = Structure.from_file(filename)
    self.structure = structure

    self.fileID = filename
    self.name = structure.composition.reduced_formula
    self.unit_cell = structure.lattice.matrix

    #species = list(structure.composition.get_el_amt_dict().keys())
    self.num_per_type = [int(x) for x in structure.composition.get_el_amt_dict().values()]
    self.species = [x.symbol for x in structure.species]
    self.num_atoms = structure.composition.num_atoms 

    self.atoms = structure.frac_coords
    #unit_conv = 180./np.pi 
    #print(np.dot(self.atoms, self.structure.lattice.inv_matrix))
    #print(structure.cartesian_coords)
    #print(arccos(np.dot(self.unit_cell[1,:],self.unit_cell[2,:].T)[0,0]/(self.cell_b*self.cell_c)))
    #print(structure.lattice.alpha)
 
    #if self.name.split()[0] == "SUPERCELL":
    #  self.is_supercell = True
    #  self.supercell_repetitions = self.name.split()[1].split('x')
    #  self.supercell_repetitions = [int(i) for i in self.supercell_repetitions]

  def save_to_poscar(self, filename):
    """
    Save the system to POSCAR using pymatgen
    """
    self.structure.to(fmt='POSCAR', filename=filename)    

  def save_to_cif(self, filename):
    """
    Save the system to cif using pymatgen
    """
    self.structure.to(fmt='cif', filename=filename)

  def save_to_xyz(self, filename):
    """
    Save the system to XYZ
    """ 
    with open( filename, 'a' ) as F:
      F = open( filename, 'a' )
      F.write( '%d\n'%self.num_atoms )
      F.write( "XYZ\n" )
      for num,row in enumerate(self.atoms):
        try:
          F.write('%s  '%self.species[num])
        except:
          F.write('X%d '%num)
        F.write( mat2str( row, "%16.10f" ) )
      F.write( "\n" )
    
  def position(self):
    """
    Return atomic positions flattened to a single long vector.
    """
    return self.atoms.reshape((1,-1))    
    
#  def supercell(self, shape):
#    """
#    Generate a supercell of the given shape
#    """
#    l,m,n = shape
#    mult = l*m*n
#
#    supercell = atomic_set()
#
#    supercell.name = "SUPERCELL %dx%dx%d %s"%(l,m,n,self.name)
#    supercell.unit_cell = multiply(self.unit_cell,array([[l,l,l],[m,m,m],[n,n,n]]))
#    supercell.recip_cell = linalg.inv(supercell.unit_cell)
#    supercell.num_atoms = mult*self.num_atoms
#    supercell.num_per_type=["%d"%(mult*int(s)) for s in self.num_per_type]
#    supercell.species=[]
#    for i,n in enumerate(self.num_per_type):
#      supercell.species.extend(n*self.species[i])
#
#    supercell.atoms = []
# 
#    for atom in array(self.atoms):
#      for displ in iterate_all_indices([l,m,n]):
#        supercell.atoms.append( atom + dot( array(displ), array(self.unit_cell) ) )
#
#    supercell.atoms = mat(supercell.atoms)
#
#    return supercell
  
  @property
  def recip_cell(self):
    return self.structure.lattice.inv_matrix

  @property
  def cell_a(self):
    return structure.lattice.a

  @property
  def cell_b(self):
    return structure.lattice.b

  @property
  def cell_c(self):
    return structure.lattice.c
  
  @property
  def cell_alpha(self, unit_conv=180./np.pi):
    return structure.lattice.alpha
    #return unit_conv*arccos(np.dot(self.unit_cell[1,:],self.unit_cell[2,:].T)[0,0]/(self.cell_b*self.cell_c))

  @property
  def cell_beta(self, unit_conv=180./np.pi):
    return structure.lattice.beta
    #return unit_conv*arccos(np.dot(self.unit_cell[0,:],self.unit_cell[2,:].T)[0,0]/(self.cell_a*self.cell_c))

  @property
  def cell_gamma(self, unit_conv=180./np.pi):
    return structure.lattice.gamma
    #return unit_conv*arccos(np.dot(self.unit_cell[0,:],self.unit_cell[1,:].T)[0,0]/(self.cell_a*self.cell_b))
  
  @property
  def fractional_atoms(self):
    return np.dot(self.atoms,self.recip_cell)
    
#  def align_to(self,reference):    
#    shift = (self.fractional_atoms - reference.fractional_atoms)
#    ishift= shift.round()
#    if abs(shift-ishift).max() > 0.4:
#      self.debug('Atomic displacement larger than .4 - possible reshuffling of atoms', LOG_WARNING)
#      print shift
#    
#    self.atoms -= dot(ishift,self.unit_cell)
