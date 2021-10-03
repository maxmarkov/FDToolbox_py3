#!/usr/bin/env python

from fdtoolbox.calculation_set import *
from fdtoolbox.utility import *
from fdtoolbox.linear_expansion import *

import sys

try:
  directory = sys.argv[1]
except:
  directory = "."
  
try:
  saxis = [float(v) for v in sys.argv[2].split('/')]  
except:
  saxis = [0., 0., 1.]

try:
  usecache = sys.argv[3] != 'nocache'
except:
  usecache = True
  
print 'Reading from directory: %s' % directory
print 'Using saxis: %s' % `saxis`
print 'Using cache: %s' % `usecache`

cs=calculation_set.read_directory(directory, saxis, usecache)

cs.set_ionic('.*/calc_.._.*[0123]')
cs.set_lattice('.*/calc_00_...')
cs.set_groundstate('.*/calc_00_0')

calculation_set.BPH_MAXDIST = 0.25
calculation_set.BPH_CORRECTION = 0.5

cs.try_fix_displacements()  
cs.try_fix_polarizations()

lin_exp = linear_expansion( cs )
lin_exp.calculate_expansion_coefficients()

me, units = lin_exp.magneto_electric_coupling(lattice=False)
print 'Clamped cell ME coupling (%s)'%units
print mat2str( me )
print 'Maximum coupling: %f %s'%(sqrt(max(linalg.eig(dot(me,me.T))[0])), units)


me, units = lin_exp.magneto_electric_coupling()
print 'Full ME coupling (%s)'%units
print mat2str( me )
print 'Maximum coupling: %f %s'%(sqrt(max(linalg.eig(dot(me,me.T))[0])), units)

evals = linalg.eig(lin_exp.B_m_n)[0]
print "Eigenvalues of force constant matrix:"
print mat2str( evals )
print 'enum\teval\tcontribution'
for i in range(lin_exp.B_m_n.shape[0]):
  inv_B_m_n = safe_inv( lin_exp.B_m_n, -3, iterate_components = [i] )[0]
  print '%d\t%f\t%s'%(i, evals[i]*cs.groundstate.volume, mat2str( (-1e4*COUPLING_TO_GAUSS*dot(lin_exp.B_m_alpha.T, dot(inv_B_m_n, lin_exp.B_m_mu) )).flatten() ) )

print lin_exp.electric_susceptibility(lattice=False), 'eps0'
print lin_exp.magnetic_susceptibility(lattice=False), '1e-8 mu0**-1' 
