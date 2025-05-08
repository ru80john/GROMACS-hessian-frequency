# GROMACS Hessian frequency
## Introduction

<div align="center">
A tool for calculating molecule vibrational mode and its corresponding frequency using GROMACS input 
and compare the results with QM calculation.
</div>

## Necessary files
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
* BEN  
```python
compound = 'BEN' # Molecule
mypath = './' # 
vib_scaling = 0.957
include_nonbonded = False
have_dihedral = False
```

![alt text](image-1.png)


some source code: https://github.com/selimsami/qforce