# Finite Differences Toolbox 

Python3 re-implementation of Finite Differences Toolbox for calculation of magnetic, dielectric, and magnetoelectric response properties of simple materials. The code is interfaced with [Pymatgen](https://pymatgen.org/#) (Python Materials Genomics) library to facilitate pre- and post-processing.

The methodology for these calculations is described in the following references:

- *"Ab Initio Indications for Giant Magnetoelectric Effects Driven by Structural Softness"* J.C. Wojdel, J. Iniguez, [PRL 105, 037208 (2010)](http://dx.doi.org/10.1103/PhysRevLett.105.037208)

- *"Magnetoelectric Response of Multiferroic BiFeO3 and Related Materials from First-Principles Calculations"* J.C. Wojdel, J. Iniguez, [PRL 103, 267205 (2009)]((http://dx.doi.org/10.1103/PhysRevLett.103.267205))

- *"First-Principles Approach to Lattice-Mediated Magnetoelectric Effects"* J. Iniguez, [PRL **101**, 117201 (2008)](http://dx.doi.org/10.1103/PhysRevLett.101.117201)

The development has been done for the following work that must be cited as well:

- *"Ferroelectricity and multiferroicity in anti–Ruddlesden–Popper structures"* M. Markov, L. Alaerts, H. Miranda, G. Petretto, W. Chen, J. George, E. Bousquet, P. Ghosez, G.-M. Rignanese, and G. Hautier, [PNAS 118 (17) e2026020118 (2021)](https://doi.org/10.1073/pnas.2026020118)

**The articles above must be cited if you use this code.**

**Note:** the code computes ion- and lattice-mediated contributions to the spin component of magneto-electric response. It does not compute the orbital component and the electronic contribution to the spin component of magneto-electric response.
