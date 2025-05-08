# %%
import math
import numpy as np
from ase.units import Hartree, mol, kJ, Bohr
import sys

class QM:
    def __init__(self, compound, mypath, vib_scaling):
        self.compound = compound
        self.mypath = mypath
        self.fchk_file = f'{self.mypath}{self.compound}/QM_GAS_PHASE/{self.compound}.fchk'
        self.vib_scaling = vib_scaling
        self.n_atoms, self.charge, self.multiplicity, self.elements, self.coords, self.hessian = self._read_fchk_file()
    
    def _read_fchk_file(self):
            n_atoms, charge, multiplicity, elements, coords, hessian = None, None, None, [], [], []
            with open(self.fchk_file, "r", encoding='utf-8') as file:
                for line in file:
                    if "Charge                                     I" in line:
                        charge = int(line.split()[2])
                    if "Multiplicity                               I" in line:
                        multiplicity = int(line.split()[2])
                    if "Number of atoms  " in line:
                        n_atoms = int(line.split()[4])
                    if "Atomic numbers  " in line:
                        n_line = math.ceil(n_atoms/6)
                        for i in range(n_line):
                            line = file.readline()
                            ids = [int(i) for i in line.split()]
                            elements.extend(ids)
                    if "Current cartesian coordinates   " in line:
                        n_line = math.ceil(3*n_atoms/5)
                        for i in range(n_line):
                            line = file.readline()
                            coords.extend(line.split())
                    if "Cartesian Force Constants  " in line:
                        n_line = math.ceil(3*n_atoms*(3*n_atoms+1)/10)
                        for i in range(n_line):
                            line = file.readline()
                            hessian.extend(line.split())
            try:
                coords = np.asfarray(coords, float)
                hessian = np.asfarray(hessian, float)
            except:
                coords = np.asarray(coords, float)
                hessian = np.asarray(hessian, float)
            coords = np.reshape(coords, (-1, 3))
            elements = np.array(elements)
            coords = coords * Bohr
            hessian = hessian * Hartree * mol / kJ / Bohr**2

            # check the format of the data
            n_atoms = self.check_type(n_atoms, 'n_atoms', int)
            charge = self.check_type(charge, 'charge', int)
            multiplicity = self.check_type(multiplicity, 'n_atoms', int)
            elements = self.check_type_and_shape(elements, 'elements', int, (n_atoms,))
            coords = self.check_type_and_shape(coords, 'coords', float, (n_atoms, 3))
            hessian = self.check_type_and_shape(hessian, 'hessian', float,
                                                 (((n_atoms*3)**2+n_atoms*3)/2,)) * self.vib_scaling**2
            
            
            return n_atoms, charge, multiplicity, elements, coords, hessian
    
    @staticmethod
    def check_type(value, name, expected_type):
        if not isinstance(value, expected_type):
            sys.exit(f'WARNING: A valid "{name}" property was not found in the hessian output'
                     ' file(s). Exiting...\n\n')
        return value

    @staticmethod
    def check_type_and_shape(value, name, expected_type, expected_shape):
        value = np.asarray(value)
        if value.size == 0:
            sys.exit(f'ERROR: No data found in the QM Hessian output file(s) for "{name}".'
                     ' Exiting...\n\n')
        elif value.dtype != np.dtype(expected_type):
            raise TypeError(f'"{name}" property expected a type of "{np.dtype(expected_type)}",'
                            f' but got "{value.dtype}" for the QM Hessian output.')
        elif value.shape != expected_shape:
            sys.exit(f'ERROR: "{name}" property expected a shape of "{expected_shape}", but got '
                     f'"{value.shape}" for the QM Hessian output. Exiting...\n\n')
        return value


