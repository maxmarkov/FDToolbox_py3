# FDToolbox_py3

Python3 re-implementation of [Finite Differences Toolbox](https://github.com/jcwojdel/FDToolbox) for calculation of magnetic, 
dielectric, and magnetoelectric response properties of simple materials. The code is interfaced with [Pymatgen](https://pymatgen.org/) 
(Python Materials Genomics) library to facilitate pre- and post-processing.

The methodology for these calculations is described in the following references:

  - Ab Initio Indications for Giant Magnetoelectric Effects Driven by
  Structural Softness
  Wojdel, Jacek C.; Iniguez, Jorge
  Physical Review Letters Vol: 105 ( 3 ) 037208 2010
  DOI: http://dx.doi.org/10.1103/PhysRevLett.105.037208

  - Magnetoelectric Response of Multiferroic BiFeO3 and Related Materials
  from First-Principles Calculations
  Wojdel, Jacek C.; Iniguez, Jorge
  Physical Review Letters Vol: 103 ( 26 ) 267205 2009
  DOI: http://dx.doi.org/10.1103/PhysRevLett.103.267205

  - First-Principles Approach to Lattice-Mediated Magnetoelectric Effects
  Iniguez, Jorge
  Physical Review Letters Vol: 101 ( 11 ) 117201  2008
  DOI: http://dx.doi.org/10.1103/PhysRevLett.101.117201
  
  These articles should be cited if you use this code. 
  
  **Note:** the code computes ion- and lattice-mediated contributions to the spin component of magneto-electric response. 
  It does not compute the orbital component and the electronic contribution to the spin component of magneto-electric 
  response. 

## Examples:

We assume that POSCAR file is in the ROOT directory.

Generate folders for the noSOC calculations: 
```
python3 examples/generate_nosoc.py
```

Compute the dielectric tensor:
```
python3 examples/compute_dielectric.py
```

Compute the magnetoelectric tensor:
```
python3 examples/compute_magnetoelectric.py
```
