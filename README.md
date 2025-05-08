# GROMACS Hessian frequency
## Introduction

<div align="center">
A tool for calculating molecular vibrational modes and corresponding frequencies using GROMACS input, and comparing the results with quantum mechanical (QM) calculations.
</div>

## Requirement
### Python Package
* numpy
* ase
* seaborn
### Necessary files
```
----FORCEFILED
|   ----molecule.ff
|       ----forcefield.itp
|       ----*.rtp
|   ----molecule.top
|
----QM_GAS_PHASE
|   ----molecule.fchk
|
----STRUCTURE
|   ----molecule.pdb
```

## Example
### BEN (benzene)  
#### Input setting
```python
compound = 'BEN' # Molecule
mypath = './' # The path can find the directory of molecule

# The scaling value of vibrational frequecy, 0.957 is used for wb97XD/6-311++G(d,p) 
# Other scaling value can be found on https://cccbdb.nist.gov/vibscalejustx.asp
vib_scaling = 0.957 

include_nonbonded = False # True: Include nonbonded interaction in calculation
have_dihedral = False # True: Include dihedral term in calculation
```
#### Results
![alt text](image-1.png)


some source code: https://github.com/selimsami/qforce