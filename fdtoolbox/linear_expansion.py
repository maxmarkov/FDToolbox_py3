import numpy as np
from fdtoolbox.calculation_set import *
from fdtoolbox.utility import invert_with_warning

class linear_expansion(): #(loggable):
  def __init__(self, calcset):
    self.symmetrize = True
    self.force_asr = False
    self.calcset = calcset
    
  def calculate_expansion_coefficients(self, update_from_calcset=True):
    if update_from_calcset is True:
      self.volume = self.calcset.groundstate.volume
      self.B_alpha_beta = np.zeros((3,3))
      self.B_alpha_mu = np.zeros((3,3))
      self.B_mu_nu = np.zeros((3,3))
       
      self.B_m_n     = self.calcset.force_constant_matrix('ion', 'ion')[0]  # force-constants matrix, climped cell
      self.B_j_k     = self.calcset.force_constant_matrix('lat', 'lat')[0]  # force-constants matrix, variable cell
      #int(self.calcset.force_constant_matrix('lat', 'ion')[0].shape)
      #print(self.calcset.force_constant_matrix('ion', 'lat')[0].T.shape)
      
      #self.B_m_j     = 0.5*( self.calcset.force_constant_matrix('lat', 'ion')[0]+ 
      #                       self.calcset.force_constant_matrix('ion', 'lat')[0].T ) # the symmetrized version of the internal strain tensor ^_{m,j} 
      # symmetrize the force constants
      if self.symmetrize:
        self.B_m_n     = 0.5*( self.B_m_n + self.B_m_n.T )
        self.B_j_k     = 0.5*( self.B_j_k + self.B_j_k.T )

      # apply the acoustic sum rule for the force-constants matrix (if activated with -a option when run)  
      if self.force_asr:
        for i in range(0,self.B_m_n.shape[0],3):
          asr = zeros([3,3])
          for j in range(0,self.B_m_n.shape[0],3):      
            asr += self.B_m_n[i:i+3,j:j+3]
      
          self.B_m_n[i:i+3,i:i+3] -= asr
        
  
  
      #self.B_m_alpha = -self.calcset.electric_polarization_matrix('ion')[0]  # the Born effective charges, climped cell
      #self.B_alpha_j = -self.calcset.electric_polarization_matrix('lat')[0]  # the Born effective charges, variable cell
      #self.B_m_mu    = self.calcset.magnetic_polarization_matrix('ion')[0]   # the magnetic charges, climped cell
      #self.B_mu_j    = self.calcset.magnetic_polarization_matrix('lat')[0]   # the magnetic charges, variable cell

    ## inverse the force-constant matrix K_{m,n}
    #inv_B_m_n = invert_with_warning(self.B_m_n, self.calcset.TRANSLATIONAL_MODE_THRESHOLD,
    #                                'Inverting force constant matrix resulted in %d soft modes', 3)      

    #self.Bhat_alpha_beta = self.B_alpha_beta -    np.dot(self.B_m_alpha.T, np.dot(inv_B_m_n, self.B_m_alpha) )   # electric_susceptibility, climped cell
    #self.Bhat_mu_nu      = self.B_mu_nu      -    np.dot(self.B_m_mu.T,    np.dot(inv_B_m_n, self.B_m_mu) )      # magnetic_susceptibility, climped cell
    #self.Bhat_j_k        = self.B_j_k        -    np.dot(self.B_m_j.T,     np.dot(inv_B_m_n, self.B_m_j) )       # elastic tensor C_jk, climped cell 
    #self.Bhat_alpha_mu   = self.B_alpha_mu   -    np.dot(self.B_m_alpha.T, np.dot(inv_B_m_n, self.B_m_mu) )      # magnetoelectric tensor, climped cell
    #self.Bhat_alpha_j    = self.B_alpha_j    -    np.dot(self.B_m_j.T,     np.dot(inv_B_m_n, self.B_m_alpha) )   # piezoelectric tensor, climped cell
    #self.Bhat_mu_j       = self.B_mu_j       -    np.dot(self.B_m_j.T,     np.dot(inv_B_m_n, self.B_m_mu) )      # piezomagnetic tensor, climped cell

    #inv_B_j_k = invert_with_warning(self.Bhat_j_k, self.calcset.ROTATIONAL_MODE_THRESHOLD,
    #                                'Inverting elastic constant matrix resulted in %d soft modes', 3)      
    #
    #self.Bgot_alpha_beta = self.Bhat_alpha_beta -    np.dot(self.Bhat_alpha_j.T, np.dot(inv_B_j_k, self.Bhat_alpha_j) ) # electric_susceptibility, variable cell
    #self.Bgot_mu_nu      = self.Bhat_mu_nu      -    np.dot(self.Bhat_mu_j.T,    np.dot(inv_B_j_k, self.Bhat_mu_j) )    # magnetic_susceptibility, variable cell
    #self.Bgot_alpha_mu   = self.Bhat_alpha_mu   -    np.dot(self.Bhat_alpha_j.T, np.dot(inv_B_j_k, self.Bhat_mu_j) )    # magnetoelectric tensor, variable cell

#  def born_charges(self):
#    return -self.volume*self.B_m_alpha, "|e|"
  
#  def magnetic_strengths(self):
#    return self.B_m_mu, "mu_B A**-1"
  
#  def piezoelectric_stress_tensor(self, ionic=True):
#    if ionic is True:
#      rval = self.Bhat_alpha_j
#    else:
#      rval = self.B_alpha_j
#    return -EA2_TO_CM2*rval, "C m**-2"
    
#  def piezoelectric_strain_tensor(self, ionic=True):
#    if ionic is True:
#      sjk = invert_with_warning(self.Bhat_j_k, self.calcset.ROTATIONAL_MODE_THRESHOLD,
#                              'Inverting elastic constant matrix resulted in %d soft modes', 3)      
#      eaj= self.Bhat_alpha_j
#    else:
#      sjk = invert_with_warning(self.B_j_k, self.calcset.ROTATIONAL_MODE_THRESHOLD,
#                              'Inverting elastic constant matrix resulted in %d soft modes', 3)      
#      eaj= self.B_alpha_j
#      
#    return -1000*EA2_TO_CM2*dot(sjk, eaj)/EVA3_TO_GPA, "pC N**-1"
  
#  def piezomagnetic_stress_tensor(self, ionic=True):
#    if ionic is True:
#      rval = self.Bhat_mu_j
#    else:
#      rval = self.B_mu_j
#      
#    return -BMA3_TO_AM*rval, "A m**-1"
    
#  def piezomagnetic_strain_tensor(self, ionic=True):
#    if ionic is True:
#      sjk = invert_with_warning(self.Bhat_j_k, self.calcset.ROTATIONAL_MODE_THRESHOLD,
#                                'Inverting elastic constant matrix resulted in %d soft modes', 3)    
#      mmj = self.Bhat_mu_j
#    else:
#      sjk = invert_with_warning(self.B_j_k, self.calcset.ROTATIONAL_MODE_THRESHOLD,
#                                'Inverting elastic constant matrix resulted in %d soft modes', 3)    
#      mmj = self.B_mu_j
#       
#    return -BMEV_TO_T*dot(sjk, mmj), "T-1"
    
#  def elastic_tensor(self, ionic=True):
#    if ionic is True:
#      rval = self.Bhat_j_k
#    else:
#      rval = self.B_j_k
#    return EVA3_TO_GPA*rval, "GPa" 
    
#  def compliance_tensor(self, ionic=True):
#    # compliance tensor is an inverted elastic tensor
#    cjk = self.elastic_tensor(ionic)[0]
#    sjk = invert_with_warning(cjk, self.calcset.ROTATIONAL_MODE_THRESHOLD*EVA3_TO_GPA,
#                              'Inverting elastic constant matrix resulted in %d soft modes', 3)
#    return 1000*sjk, "TPa**-1"  
  
#  def magneto_electric_coupling(self, lattice=True):
#    if lattice is True:
#      rval = self.Bgot_alpha_mu
#    else:
#      rval = self.Bhat_alpha_mu
#
#    return -1e4*COUPLING_TO_GAUSS*rval, "1e-4 g.u."
  
#  def force_constant_matrix(self):
#    return self.calcset.groundstate.volume*self.B_m_n, "eV Andstrom**-2  "
  
  def electric_susceptibility(self, lattice=True):
    if lattice is True:
      rval = self.Bgot_alpha_beta
    else:
      rval = self.Bhat_alpha_beta
    return(-EEAEV_TO_EPSILON0*rval, "epsilon0")
      
#  def magnetic_susceptibility(self, lattice=True):
#    if lattice is True:
#      rval = self.Bgot_mu_nu
#    else:
#      rval = self.Bhat_mu_nu
#    return -1e8*BM2A3EV_TO_1MU0*rval, "1e-8 mu0**-1"
