import numpy as np
from pymatgen.io.vasp.outputs import Outcar, Vasprun
from pymatgen.core.periodic_table import Element

from fdtoolbox.utility import *
#from fdtoolbox.atomic.atomic_set import atomic_set
#from fdtoolbox.atomic.utility import loggable

import os
import atexit
import re
import sys
try:
  import gzip
  have_gzip = True
except:
  have_gzip = False

class calculation():
  """
  A class for keeping all data read from OUTCAR for a single calculation. It reads the following
  properties of the calculation and stores them in the following attributes:
  - self.unit_cell - unit cell (3x3 A)
  - self.recip_cell - reciprocal unit cell (3x3 A**-1)
  - self.volume - volume calculated from the unit cell (1 A**3)
  - self.energy - total energy reported by VASP (1 eV)
  - self.name - name of the system as put in POSCAR
  - self.num_per_type - list of species numbers 
  - self.num_atoms - number of atoms
  - self.atoms - atomic positions internal coordinates (self.num_atomsx3 1)
  - self.forces - forces acting on ions as energy per internal coordinate derivative (self.num_atomsx3 eV)
  - self.stress - stress on the cell in cartesian coordinates (3x3 eV/A**3)
  - self.magnetization - total cell magnetization in cartesian coordinates (3 e)
  - self.polarization - polarization as read by load_polarization() method (3 e/A**2)
  For other polarization related attributes see load_polarization() method documentation.
  """

  saxis = [ 0., 0., 1. ]
  
  def load_polarization(self, filename):
    """
    This function tries to read Berry phase calculations from three files:
    filename+'_berry_1' filename+'_berry_2' filename+'_berry_3'
    If the files are not readable, no polarization attribute is set.
    On successful run this method fills the following attributes:
    - self.polarization - polarization as calculated from the three files, in cartesian coordinates (3 e/A**2)
    - self.ev_phase - expectation value reported by VASP for each of the three directions (3x3 -e*A)
    - self.ev_recip - the mean of the expectation value projected to internal coordinates (3 -e)
    - self.berry_phase - berry phases of all strings in each of the directions (3xN 1)
    - self.polarization_quant - polarization quant for the calculation (3x3 e/A**2)
    - self.i_polarization_quant - inverse of polarization quant
    """
    filename = filename[:-7]
    print(filename)
    self.berry_phase=3*[None]
    self.berry_ev=3*[None]
    self.berry_ion=3*[None]
    self.num_per_type = []
    for kdirection in range(1,4):
      fname = os.path.join(filename, 'Berry_%d/OUTCAR' % kdirection)
      if os.access( fname, os.R_OK ):
        ifile = open(fname)
  #    else:
  #      # Make sure, we have 0,0,0 as polarisation in case there is no berry phase output available 
  #      self.berry_phase=3*[array([[0.,0.]])]
  #      self.berry_ev=3*[array([0., 0., 0.])]
  #      self.berry_ion=3*[[[0., 0., 0.]]]
  #      self.polarization_quant = self.unit_cell.copy() / self.volume
  #      self.berry_multiplier = 1.
  #      self.polarization_quant *= self.berry_multiplier
  #      self.i_polarization_quant = linalg.inv(self.polarization_quant)
  #      self.ev_recip = dot( mean(self.berry_ev,axis=0), self.recip_cell )
  #      self.recalculate_polarization()
  #      return      
      line = " "
      self.berry_phase[kdirection-1] = []
      self.berry_ev[kdirection-1] = []
      while line:
        line = ifile.readline()
        if line.startswith("   ISPIN  ="):
          data = line.split()
          self.numspin = int(data[2])
        elif line.startswith("   LNONCOLLINEAR ="):
          # This is a bit tricky - noncollinear calculations are officially numspin==1, but
          # they do contain separate electrons, not electron pairs, and therefore should be 
          # treated as numspin==2 for calculation of polarizations.
          data = line.split()
          if data[2] == "T":
            self.numspin = 2
        elif line.startswith("Expectation value"): #VASP 4.x
          line = ifile.readline()
          data = line.split() 
          self.berry_ev[kdirection-1].append(  [ float(data[3][:-1]), float(data[4][:-1]), float(data[5]) ] )
        elif line.find("e<r>_ev") > 0: #VASP 5.x
          data = line.split()
          self.berry_ev[kdirection-1].append( [ float(data[-5]), float(data[-4]), float(data[-3]) ] )
        elif line.startswith("K-point string"): #VASP 4.x 5.x
          data = line.split()
          last_kpoint_weight = float(data[5])
        elif line.startswith("Berry-Phase term:  "): #VASP 4.x
          data = line.split()
          self.berry_phase[kdirection-1].append( [float(data[2]), last_kpoint_weight] )
        elif line.startswith("              Im ln[Det|M_k|]=" ): #VASP 5.x
          data = line.split()
          self.berry_phase[kdirection-1].append( [float(data[-2]), last_kpoint_weight] )
        #No need to read ionic term - we will recalculate it anyway   
        #elif line.startswith("ionic term:"): #VASP 4.x
        #  line = ifile.readline()
        #  data = line.split() 
        #  self.berry_ion[kdirection-1] = [ [ float(data[3][:-1]), float(data[4][:-1]), float(data[5]) ] ]
     
      print('numspin =', self.numspin)
      print('berry_ev', self.berry_ev)
      print('berry_phase', self.berry_phase)   
      self.berry_phase[kdirection-1] = array(self.berry_phase[kdirection-1])
      self.berry_ev[kdirection-1] = mean(array(self.berry_ev[kdirection-1]), axis=0)
          
    self.polarization_quant = self.unit_cell.copy() / self.volume
            
    # In case of nonpolarised calculations, we have two electrons per each Berry string
    # So the quantum needs to be doubled.
    self.berry_multiplier = 3-self.numspin
    self.polarization_quant *= self.berry_multiplier
    self.i_polarization_quant = np.linalg.inv(self.polarization_quant)
    
  #  self.ev_recip = np.dot( np.mean(self.berry_ev,axis=0), self.recip_cell )
  #  self.recalculate_polarization()
    
  #def berry_term(self):
  #  """
  #  Returns Berry term part of the polarization for each of the IGPAR directions (i.e. 
  #  k-point weighted sum of Berry terms for each k-point). The term is multiplied by proper 
  #  multiplier to make sure that closed-shell calculations sum up to proper number of electrons. 
  #  """
  #  bp = zeros((1,3))
  #  for i in range(3):
  #    # berry_phase[i][:,0] -> berry term, self.berry_phase[i][:,1] -> k-point weight
  #    bp[0,i] = sum(self.berry_phase[i][:,0]*self.berry_phase[i][:,1])
  #  return self.berry_multiplier*mat(bp)
    
  #def ionic_term(self):
  #  """
  #  Returns ionic part of polarization. Contrary to VASP it calculates this with 0 as a 
  #  reference point (VASP uses fract(1/2,1/2,1/2)).
  #  """
  #  ion=zeros((1,3))
  #  for atom, zval in zip(self.atoms, self.zvals):
  #    ion += -zval*atom
  #  return mat(ion)
    
      
  #def recalculate_polarization(self):
  #  """
  #  Recalculates the polarization from the data gathered from the OUTCAR_berry_[123] files.
  #  This function should be called whenever we tamper with Berry phases (to fix problems with
  #  polarization quantum).
  #  """    
  #  # Gather data from all three directions

  #  self.polarization = self.berry_term()*self.unit_cell + mean(self.ionic_term(),axis=0) + mean(self.berry_ev,axis=0)

  #  # Note that VASP reports data in 'electrons * A', meaning the lack of minus sign
  #  # and not scaling by cell volume.
  #  self.polarization /= -self.volume

  #  return self.polarization
  #  # At the end, the polarization is expressed in e/A**2
    


  def load_from_outcar(self, filename, structure):
    """
    This function reads atom positions, unit cell, forces and stress tensor from the VASP
    OUTCAR file
    """
    if os.access(filename, os.R_OK):
      outcar = Outcar(filename)
    else:
      raise Exception("Missing OUTCAR file in {}".format(filaname))
      
    self.name = 'rr'
 
    # ATOMIC PROPERTIES
    self.num_atoms = int(structure.composition.num_atoms)                                     # total number of atoms
    self.atoms_frac = structure.frac_coords                                                   # fraction coordinates
    self.atoms_cart = structure.cart_coords                                                   # cartesian coordinates
    self.num_per_type = [int(x) for x in structure.composition.get_el_amt_dict().values()]    # number of atoms per type
    self.species = [x.symbol for x in structure.species]                                      # list with atomic symbols for each specie
    self.pomass = [Element(x.symbol).atomic_mass for x in structure.species]                  # list with atomic masses for each specie
    self.charges = [x['tot'] for x in outcar.charge]                                          # list with charges
    
    # read zval dict and create a mapping 
    outcar.read_pseudo_zval()
    zval_dict = outcar.zval_dict
    self.zvals = [zval_dict[x.symbol] for x in structure.species]
   
    # UNIT CELL
    self.unit_cell = structure.lattice.matrix
    self.volume = structure.lattice.volume

    # total energy in eV
    self.energy = outcar.final_energy
    # atomic positions read from Outcar
    self.atoms = np.array(outcar.read_table_pattern(header_pattern=r"\sPOSITION\s+TOTAL-FORCE \(eV/Angst\)\n\s-+",
                                                     row_pattern=r"\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+",
                                                     footer_pattern=r"\s--+",
                                                     postprocess=lambda x: float(x),
                                                     last_one_only=False)[0])
    # forces acting on atoms from Outcar
    self.forces = np.array(outcar.read_table_pattern(header_pattern=r"\sPOSITION\s+TOTAL-FORCE \(eV/Angst\)\n\s-+",
                                                     row_pattern=r"\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)",
                                                     footer_pattern=r"\s--+",
                                                     postprocess=lambda x: float(x),
                                                     last_one_only=False)[0])

    # READ STRESS: (XX, YY, ZZ, XY, YZ, ZX) and convert values to array
    outcar.read_pattern({'stress':  r"in kB\s+([\.\-\d]+)\s+([\.\-\d]+)\s+([\.\-\d]+)\s+([\.\-\d]+)\s+([\.\-\d]+)\s+([\.\-\d]+)"}, terminate_on_match=False, postprocess=float)
    stress = outcar.data.get("stress")[0]
    self.stress = np.array([[stress[0], stress[3], stress[5]],
                            [stress[3], stress[1], stress[4]],
                            [stress[5], stress[4], stress[2]]])

    # MAGNETIZATION and its PROJECTION ON EVERY ATOM
    try:
       # non-collinear
       outcar.read_pattern({'total_mag': r"number of electron\s+\S+\s+magnetization\s+([\.\-\d]+)\s+([\.\-\d]+)\s+([\.\-\d]+)\s"}, terminate_on_match=False, postprocess=float)
       self.__magnetization = outcar.data.get("total_mag")[-1]

       self.proj_magn = np.array([mag['tot'].moment for mag in outcar.magnetization]) # 2-D array
    except:
       # collinear
       outcar.read_pattern({'total_mag': r"number of electron\s+\S+\s+magnetization\s+(" r"\S+)"}, terminate_on_match=False, postprocess=float)
       self.__magnetization = np.zeros(3)
       self.__magnetization[2] = outcar.data.get("total_mag")[-1][0]
 
       self.proj_magn = np.zeros((self.num_atoms,3))
       for i in range(self.num_atoms):
           self.proj_magn[i,2] = outcar.magnetization[i]['tot']                       # 2-D array

    # NOTE: sum of magnetization projected on atoms differs from the total magnetization
    self.proj_magn_sum = np.sum(self.proj_magn, axis=0)

    # POLARIZATION    
    self.load_polarization(filename)
    self.fileID = filename

    
  #def get__magnetization(self):
  #  """
  #  Returns the magnetization vector in Cartesian coordinates. The self.__magnetization always
  #  stores the magnetization as read from OUTCAR. Note that for optimization purposes the rotated
  #  magnetization is stored between calls. If self.saxis is changed, the invalidate_magnetization 
  #  method _must_ be called on all Calculation class objects.    
  #  """
  #  # Rotation of reported magnetic moment - according to:
  #  # http://cms.mpi.univie.ac.at/vasp/vasp/node159.html
  #  if not hasattr(self, "__rotated_magnetization"):
  #    self.__rotated_magnetization = dot( self.__magnetization, inv_saxis_rotation(self.saxis) )
  #    
  #  return self.__rotated_magnetization 
  
  #def invalidate_magnetization(self):
  #  """
  #  Clear magnetization expressed in cartesian coordinates.
  #  """
  #  del self.__rotated_magnetization
  #  
  #magnetization = property( get__magnetization, None )
    
  def force(self, mode = 'ion'):
    """
    Return forces acting on 'ion's, 'lat'tice, or 'all' of them, flattened to one long vector.
    Note that the forces on lattice are expressed in terms of change of total energy of the 
    system per (cartesian) strain - which brings them to the same units as forces acting on ions (eV/1).
    """
    if mode == 'ion':
      return(self.forces.reshape((1,-1)) / self.volume)
    elif mode == 'lat':
      return(self.stress.flatten(0))
    else:
      ion = self.force('ion')
      lat = self.force('lat')
      return(bmat('ion, lat'))

  def projected_charges(self):
    if hasattr(self, "charges"):
      return(array(self.zvals)-array(self.charges))
    else:
      return(None)
  
  #def projected_magnetizations(self):
  #  if hasattr(self, "proj_magn"):
  #    return dot( self.proj_magn, inv_saxis_rotation(self.saxis) )
  #  else:
  #    return None
    
  
class calculation_set():
  """
  A class representing a full set of calculations read from multiple subdirectories.
  """

  """
  These two class attributes define the modulo factor by which the berry phases are regulated when 
  aligning the phases of different calculations. Normally they should be left alone. However, in case of 
  mixed set of calculations (ISPIN=1 and ISPIN=2), the way the different spins are treated might lead to 
  a shift by 0,5 of the berry phases reported for one of the spins. The standard way of aligning to the 
  first phase of groundstate calculation will likely fail (as beta-spin electrons will be contributing 
  at exactly the noncontinguous +/-0.5 region). A quick fix to this is to divide both values by half.
  A better one would be to remember the different spin components...
  """
  BPH_MAXDIST = 0.5
  BPH_CORRECTION = 1.0

  
  # Used as a threshold for seeking if the positions of polarizations need to be fixed.
  MAX_DISPLACEMENT = 0.1
  MAX_POLARIZATION_SHIFT = 0.1

  # Threshold used to recognize translational modes of force constant matrix
  TRANSLATIONAL_MODE_THRESHOLD = 0.001
  
  # Threshold used to recognize rotational modes of elastic constant matrix
  ROTATIONAL_MODE_THRESHOLD = 0.001

  # The stiffness of the constraint on lattice distortions used in add_lattice_constraint()
  # and ionic distortions in me_coupling. 
  CONSTRAINT_STIFFNESS = 1e8
  
  # Mapping between descriptive names and lists of calculations read.
  mode_mapping={'all': '_calculation_list',
                'ion': '_ionic_displacements_list',
                'lat': '_lattice_displacements_list'}
  
  def __init__(self):
    self._calculation_list=[]
    self._ionic_displacements_list=[]
    self._lattice_displacements_list=[]
  
  @classmethod
  def read_directory(cls, dirname, fdt, saxis, usecache=False):
    """
    Read the whole directory containing subdirectories with OUTCARS.
    This function fills in the self._calculation_list attribute.
    At the end the self.groundstate attribute is chosen from the list 
    based on the energy.
    """
    import pickle
  
  #calculation.saxis = saxis
   
  #  if usecache:
  #    dataname = dirname+'_pickled.dat'
  #    if os.access( dataname, os.R_OK ):
  #      F = open(dataname,'rb')
  #      cs = pickle.load(F)
  #      F.close()
  #      return cs
    
    cs = calculation_set()

    for f in fdt.structures_ion:
      path = os.path.join(dirname,f)
      if os.path.isdir(path):
        fpath = os.path.join(path,'OUTCAR')
        if os.access( fpath, os.R_OK ) or os.access( fpath+'.gz', os.R_OK ):
          print(fpath)
          c = calculation()
  #        cs.debug('Reading calculation from %s' % fpath, LOG_INFO)
          c.load_from_outcar(fpath, fdt.structures_ion[f])
  #        cs._calculation_list.append( c )
      break

  #  me = 0.0
  #  for c in cs._calculation_list:
  #    if c.energy < me:
  #      me = c.energy
  #      cs.groundstate = c
      
  #  supercell_definition = re.search("(SUPERCELL)\s+([^\s]+)", cs.groundstate.name) 
  #  if supercell_definition:
  #    print supercell_definition.groups()
  #    repetitions = re.match("([0-9]+)x([0-9]+)x([0-9]+)", supercell_definition.groups()[1])
  #    if repetitions:
  #      cs.expand_translational_symmetry([int(r) for r in repetitions.groups()])
  #    else:
  #      #TODO: load a file with translationlist?
  #      cs.expand_translational_symmetry(xxx)
        
  #  if usecache:
  #    F = open(dataname,'wb')
  #    pickle.dump(cs,F,-1)
  #    F.close()

  #  return cs

#  #def add_calculation(self, c ):
#  #  self._calculation_list.append( c )
#    
#  #def try_fix_displacements(self):
#  #  """
#  #  Fix displacements within the calculation set. Ions are moved around by 
#  #  lattice displacements in order to fit within self.MAX_DISPLACEMENT 
#  #  from the groundstate (so the groundstate must be properly set before 
#  #  calling this function).
#  #  """
#  #            
#  #  def fix_one_calculation( item ):
#  #    for idx,atom in enumerate(item.atoms):
#  #      gs_atom = self.groundstate.atoms[idx]
#  #      for dx,dy,dz in iterate_shifts():
#  #        new_atom = atom+[[dx,dy,dz]]*item.unit_cell
#  #        if metro_dist(new_atom, gs_atom) <= self.MAX_DISPLACEMENT:
#  #          self.debug('    moving atom %d by %d,%d,%d' % (idx,dx,dy,dz),LOG_ALLINFO)
#  #          item.atoms[idx] = new_atom[0]
#  #    
#
#  #  gs_position = self.groundstate.position()
#  #  for item in self.calculations('ion'):
#  #    if metro_dist(item.position(), gs_position) > self.MAX_DISPLACEMENT:
#  #      self.debug("   fixing %s" % item.fileID, LOG_ALLINFO)
#  #      fix_one_calculation( item )
#  #      item.recalculate_polarization()
#  #        
#  #  self.debug("Corrected displacements", LOG_ALLINFO)
#  #  self.debug(mat2str(self.displacements('ion')), LOG_ALLINFO)
#
#
#  #def try_fix_polarizations(self):
#  #  """
#  #  Try to fix the polarization if it is screwed by polarization quantum.
#  #  """
#  #  if self.groundstate.berry_phase[0] is None:
#  #    return
#  #  
#  #  def fix_berry_phases( calc ):
#  #    changed = False
#  #    for kdir in range(3):
#  #      reference_phase = self.groundstate.berry_phase[kdir][0,0]
#  #      for ph in calc.berry_phase[kdir]:
#  #        while ph[0]-reference_phase > calculation_set.BPH_MAXDIST:
#  #          ph[0] -= calculation_set.BPH_CORRECTION
#  #          changed = True
#  #        while reference_phase-ph[0] > calculation_set.BPH_MAXDIST:
#  #          ph[0] += calculation_set.BPH_CORRECTION
#  #          changed = True
#  #    return changed
#  #      
#  #  for calc in self._calculation_list:
#  #    #For some reason this does not work as expected and it's safer to just fix/recalculate all calculations
#  #    #if fix_berry_phases(calc):
#  #      #p = calc.recalculate_polarization()
#  #      #self.debug('Corrected polarization %s = %f %f %f'%(calc.fileID,p[0,0],p[0,1],p[0,2]), LOG_ALLINFO)
#  #    fix_berry_phases(calc)
#  #    p = calc.recalculate_polarization()
#  #    self.debug('Corrected polarization %s = %f %f %f'%(calc.fileID,p[0,0],p[0,1],p[0,2]), LOG_ALLINFO)
#    
#  #def try_polarizations_chain(self):
#  #  def align_berry_phases( ref_calc, dest_calc ):
#  #    for kdir in range(3):
#  #      if ref_calc.berry_phase[kdir].shape != dest_calc.berry_phase[kdir].shape:
#  #        self.debug("Cannot align berry phases", LOG_WARNING)
#  #      else:
#  #        #print "IN"
#  #        #print ref_calc.berry_phase[kdir]
#  #        #print dest_calc.berry_phase[kdir]
#  #        count = 0
#  #        for i in range(ref_calc.berry_phase[kdir].shape[0]):
#  #          while ref_calc.berry_phase[kdir][i,0] - dest_calc.berry_phase[kdir][i,0] > self.BPH_MAXDIST:
#  #            dest_calc.berry_phase[kdir][i,0] += self.BPH_CORRECTION
#  #            count+=1
#  #          while ref_calc.berry_phase[kdir][i,0] - dest_calc.berry_phase[kdir][i,0] < -self.BPH_MAXDIST:
#  #            dest_calc.berry_phase[kdir][i,0] -= self.BPH_CORRECTION
#  #            count+=1
#  #        if count > 0:
#  #          print 'Number of flips: %d'%count
#  #        #print "OUT"
#  #        #print dest_calc.berry_phase[kdir]
#  #        
#  #  def align_ions( ref_calc, dest_calc ):
#  #    for idx,atom in enumerate(dest_calc.atoms):
#  #      gs_atom = ref_calc.atoms[idx]
#  #      for dx,dy,dz in iterate_shifts():
#  #        new_atom = atom+[[dx,dy,dz]]*dest_calc.unit_cell
#  #        if metro_dist(new_atom, gs_atom) <= 0.5:
#  #          loggable.debug('    moving atom %d by %d,%d,%d' % (idx,dx,dy,dz),LOG_ALLINFO)
#  #          dest_calc.atoms[idx] = new_atom[0]
#  #  
#  #  self._calculation_list[0].recalculate_polarization()
#  #  for i in range(1,len(self._calculation_list)):
#  #    print "align", self._calculation_list[i-1].fileID, self._calculation_list[i].fileID
#  #    align_berry_phases(self._calculation_list[i-1],self._calculation_list[i])
#  #    align_ions(self._calculation_list[i-1],self._calculation_list[i])
#  #    self._calculation_list[i].recalculate_polarization()
#    
#
#  #def set_groundstate(self, pattern):
#  #  """
#  #  Set the self.groundstate attribute to a given file pattern. This is usefull for constant 
#  #  volume calculations, where calculated energy is not be minimal with respect to lattice deformations 
#  #  at the optimised structure.
#  #  """
#  #  for calc in self._calculation_list:
#  #    if re.match(pattern,calc.fileID):
#  #      self.groundstate = calc
#  #      return
#  #    
#  #  self.debug("Cannot set groundstate to %s" % pattern, LOG_WARNING)
#      
#  #def calculations(self,mode='all'):
#  #  """
#  #  Returns only a subset of calculations. In this way one can easily work on 'ion'ic or 'lat'tice
#  #  distortions and treat them separately where neccessary.
#  #  """
#  #  if mode == 'ion':
#  #    return self._ionic_displacements_list
#  #  elif mode == 'lat':
#  #    return self._lattice_displacements_list
#  #  else:
#  #    return concatenate( [ self._ionic_displacements_list, self._lattice_displacements_list ] )
#  
#  #def set_subset(self,mode,pattern):
#  #  """
#  #  Set the subset of calculations based on the file pattern.
#  #  """
#  #  if mode == 'all':
#  #    return
#  #  
#  #  current_attr = self.calculations(mode)
#  #  del current_attr[:]
#  #  for c in self._calculation_list:
#  #    if re.match(pattern,c.fileID):
#  #      current_attr.append(c)
#    
#  #def set_ionic(self, pattern):
#  #  self.set_subset('ion',pattern)
#    
#  #def set_lattice(self, pattern):
#  #  self.set_subset('lat',pattern)
#      
#  #def displacements(self,mode='ion', calcset='ion'):
#  #  """
#  #  Returns displacements from the groundstate for all the calculations from the list.
#  #  This method returns fractional coordinate displacements of ions in case of mode='ion' 
#  #  and strain in case of mode='lat' for a specified calculation subset. In case of mode='all' 
#  #  it returns full generalised displacement matrix.
#  #  """
#  #  if mode == 'ion':
#  #    return mat(concatenate([ calc.position() - self.groundstate.position() for calc in self.calculations(calcset) ]))
#  #  elif mode == 'lat':
#  #    return mat(concatenate([ (self.groundstate.recip_cell*calc.unit_cell - eye(3)).flatten(1) for calc in self.calculations(calcset) ]))
#  #  else:
#  #    return hstack([self.displacements('ion', calcset), self.displacements('lat', calcset)])
#    
#  #def forces(self, mode='ion', calcset='ion'):
#  #  """
#  #  Returns forces exerted on 'mode' subset of coordinates, as read from 'calcset' subset of all calculations. 
#  #  In new formalism, all global properties are defined as per unit volume and that stresses are already represented like this.
#  #  """
#  #  return mat(concatenate([ calc.force(mode) for calc in self.calculations(calcset)]))
#  
#  #def magnetizations(self,calcset='all'):
#  #  """
#  #  Returns magnetizations read from 'calcset' subset of all calculations. 
#  #  In new formalism, all global properties are defined as per unit volume
#  #  """
#  #  return mat(concatenate([ c.magnetization for c in self.calculations(calcset)])) / self.groundstate.volume
#  
#  #def magnetization_differences(self,calcset='all'):
#  #  """
#  #  Returns magnetization differences between calcset of calculation and the groundstate.
#  #  """
#  #  #if calcset == 'ion':
#  #  return mat(concatenate([ c.magnetization-self.groundstate.magnetization for c in self.calculations(calcset)]))
#  #  #elif calcset == 'lat':
#  #  #  return mat(concatenate([ c.magnetization*c.recip_cell*self.groundstate.unit_cell-self.groundstate.magnetization for c in self.calculations(calcset)]))
#  
#  #def polarizations(self,calcset='all'):
#  #  """
#  #  Returns polarizations read from 'calcset' subset of all calculations. 
#  #  """
#  #  return concatenate([ c.polarization for c in self.calculations(calcset)])
#  
#  #def polarization_differences(self,calcset='all'):
#  #  """
#  #  Returns polarization differences between calcset of calculation and the groundstate.
#  #  Note that in case of lattice distortions, comparing the polarizations directly does not 
#  #  make much sense as they differ by the factor induced by volume cell differences and cell
#  #  rotations. Therefore this function returns difference of berry-phases multiplied by 
#  #  groundstate unit cell for lattice distortion subset of calculations.
#  #  This implements the derivative of berry phases over strain formula (24) from 
#  #  D. Vanderbilt, J Phys Chem Solid 61 (2000) 147-151   
#  #  """
#  #  if calcset=='ion':
#  #    return concatenate( [ calc.polarization-self.groundstate.polarization for calc in self.calculations(calcset) ] )
#  #  elif calcset=='lat':
#  #    return concatenate( [ -(calc.berry_term()+calc.ev_recip-self.groundstate.berry_term()-self.groundstate.ev_recip)*self.groundstate.unit_cell/self.groundstate.volume for calc in self.calculations('lat') ] )
#  #  else:
#  #    return vstack([self.polarization_differences('ion'),self.polarization_differences('lat')])
#  
#  #def force_constant_matrix(self, mode='ion', calcset='ion'):
#  #  """
#  #  Calculate force constant matrix (A) and residual forces (F0) at groundstate from solving set of linear
#  #  equations involving forces (F) and displacements (d):
#  #  F = d*A + F0
#  #  The set of equations might be overdefined - the found solution is optimal in the leastsquare sense.
#  #  
#  #  Second returned value is the residual of the force at groundstate (which should be 0)
#  #  """
#  #  
#  #  # FIXME: this is a workaround for quite strange behaviour...
#  #  # When given all the calculations, the lattice part of the force constant
#  #  # matrix becomes asymmetric and generally noisy. This does not happen when
#  #  # only subsets are given. Even off-diagonal (ion-lattice) coefficients behave
#  #  # very well, when only calculated on separate calculation sets.
#  #  if calcset == 'all':
#  #    return vstack( [ self.force_constant_matrix(mode, 'ion')[0],
#  #                     self.force_constant_matrix(mode, 'lat')[0] ]), 0
#  #  
#  #  displacements = self.displacements(calcset, calcset)
#  #  forces = self.forces(mode, calcset)
#
#  #  # min || displacemnent*x - forces||^2 
#  #  # Returns: fc_matrix - the least-square solution (force constant matrix)
#  #  #                 rs - sums of residuals: || displacemnent*fc_matrix - forces||^2    
#  #  #                 rk - rank of displacement matrix
#  #  #                  s - singular values of displacement matrix 
#  #  (fc_matrix, rs, rk, s) = linalg.lstsq( 
#  #    hstack( [ displacements, ones((displacements.shape[0],1)) ] ), forces ) 
#  #  residual = fc_matrix[-1,:]   # set residual to be the last column of force constant matrix
#
#  #  # remember that when we talk about force constant matrix, we think actually second derivative
#  #  # of the energy... which is exactly opposite sign to derivative of forces 
#  #  fc_matrix = -fc_matrix[:-1,:]
#  #  if hasattr(self,'_lattice_constraint') and (mode=='all' or calcset=='all' or (mode=='lat' and calcset=='lat')):
#  #    fc_matrix[-9:,-9:] += self.CONSTRAINT_STIFFNESS*self._lattice_constraint 
#  #  return fc_matrix, residual
#         
#  #def magnetic_polarization_matrix(self, mode='ion'):
#  #  """
#  #  Returns numerically computed derivative of total magnetization vector over general coordinates.
#  #  """
#  #  if mode == 'ion' or mode == 'lat':
#  #    displacements = self.displacements(mode, mode)
#  #    mag_dif = self.magnetization_differences(mode)
#  #    (mag_pol_mat, rs, rk, s) = linalg.lstsq( 
#  #                               hstack( [ displacements, ones((displacements.shape[0],1)) ] ), mag_dif )
#  #    residual = mag_pol_mat[-1,:]
#  #    #mag_pol_mat = mag_pol_mat[:-1.:] MM
#  #    mag_pol_mat = mag_pol_mat[:-1:]
#  #    return mag_pol_mat/self.groundstate.volume, residual
#  #  else:
#  #    impl,ir = self.magnetic_polarization_matrix('ion')
#  #    lmpl,lr = self.magnetic_polarization_matrix('lat')
#  #    return vstack([impl, lmpl]), ir+lr
#
#  #def electric_polarization_matrix(self, mode='ion'):
#  #  """
#  #  Returns numerically computed derivative of polarization vector over general coordinates.
#  #  """
#  #  if mode == 'ion' or mode == 'lat':
#  #    displacements = self.displacements(mode, mode)
#  #    pol_dif = self.polarization_differences(mode)
#  #    # min || displacements*x - dP||^2
#  #    (el_pol_mat, rs, rk, s) = linalg.lstsq( 
#  #                              hstack( [ displacements, ones((displacements.shape[0],1)) ] ), pol_dif )
#  #    residual = el_pol_mat[-1,:]
#  #    el_pol_mat = el_pol_mat[:-1,:]
#  #  
#  #    if (mode == 'ion'):
#  #      return el_pol_mat, residual
#  #    else:
#  #      return el_pol_mat, residual
#  #  else:
#  #    iepl,ir = self.electric_polarization_matrix('ion')
#  #    lepl,lr = self.electric_polarization_matrix('lat')
#  #    return vstack([iepl, lepl]), ir+lr
#  #
#  #def clear_lattice_constraints(self):
#  #  """
#  #  Clear all imposed lattice distortion constraints.
#  #  """
#  #  if hasattr(self,'_lattice_constraint'):
#  #    del self._lattice_constraint
#  #    
#  #def add_lattice_constraint(self, dir_plane_mat, flags=['xyshear','oopshear']):
#  #  """
#  #  Constrain the lattice from deforming in a given cartesian direction.
#  #  """
#  #  direction = mat(dir_plane_mat)
#  #  
#  #  if direction.shape[0] == 1:
#  #    # We have a vector, which defines uncompressible direction
#  #    direction /= sqrt(direction*direction.T)
#  #    strain = direction.T*direction
#  #  elif direction.shape[0] == 2:
#  #    # We have two vectors which define rigid plane
#  #    direction[0] /= sqrt(direction[0]*direction[0].T)
#  #    direction[1] /= sqrt(direction[1]*direction[1].T)
#  #    self.add_lattice_constraint(direction[0])                # XX
#  #    self.add_lattice_constraint(direction[1])                # YY
#  #    if 'xyshear' in flags:
#  #      self.add_lattice_constraint(direction[1].T*direction[0]) # YX
#  #      self.add_lattice_constraint(direction[0].T*direction[1]) # XY
#  #    out_of_plane = cross(direction[0],direction[1])
#  #    self.add_lattice_constraint(direction[0].T*out_of_plane) # XZ
#  #    self.add_lattice_constraint(direction[1].T*out_of_plane) # YZ
#  #    if 'oopshear' in flags:
#  #      self.add_lattice_constraint(out_of_plane.T*direction[0]) # ZX
#  #      self.add_lattice_constraint(out_of_plane.T*direction[1]) # ZY
#  #    return
#  #  elif direction.shape[0] == 3:
#  #    # We have directly strain component generalised direction
#  #    strain = direction
#  #    
#  #  strain = mat(strain.flatten(0))
#  #  constraint = strain.T*strain
#  #  
#  #  if not hasattr(self,'_lattice_constraint'):
#  #    self._lattice_constraint = mat(zeros((9,9)))
#  #    
#  #  self._lattice_constraint += constraint
#  #  
#  #def expand_translational_symmetry(self, supercell):
#  #  #FIXME!!!
#  #  # more bad things than I can really tell....
#  #  print self.groundstate.unit_cell
#  #  original_cell = self.groundstate.unit_cell / array([supercell, supercell, supercell]).T
#  #  print original_cell
#  #  
#  #  def shift_calc(calc, shift):
#  #    newcalc = copy.deepcopy(calc)
#  #    for i in range(calc.atoms.shape[0]):
#  #      newcalc.atoms[i,:] = calc.atoms[i,:]+shift
#  #    return newcalc
#
#  #  newcalclist = []
#  #  for i,j,k in iterate_all_indices(supercell):
#  #    shift = dot([i,j,k], original_cell)
#  #    for calc in self._calculation_list:
#  #      newcalclist.append(shift_calc(calc, shift))
#  #  
#  #  self._calculation_list = newcalclist
#  #  print mat2str(vstack([c.atoms for c in self._calculation_list]))
#    
