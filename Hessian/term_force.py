from Hessian.forces import calc_bonds, calc_angles, calc_rb_diheds, calc_imp_diheds, calc_pairs

class VDW_Q_F:
    def __init__(self, params):
        self.params = params

    def calc_forces(self, coords, force):
        for pair, params in self.params.items():
            calc_pairs(coords, pair, params, force[0])

class Bond_F:
    def __init__(self, params, atomlist):
        self.params = params
        self.atomlist = atomlist

    def calc_forces(self, coords ,force):
        for bond, [fconst, b_length] in self.params.items():
            atoms = bond.split('-')
            atomidx = [self.atomlist.index(atoms[0]), self.atomlist.index(atoms[1])]
            calc_bonds(coords, atomidx, b_length, fconst, force[0])

class Angle_F:
    def __init__(self, params, atomlist):
        self.params = params
        self.atomlist = atomlist  
      
    def calc_forces(self, coords, force):
        for angle, [fconst, ang] in self.params.items():
            atoms = angle.split('-')
            atomidx = [self.atomlist.index(atoms[0]), self.atomlist.index(atoms[1]), self.atomlist.index(atoms[2])]
            calc_angles(coords, atomidx, ang, fconst, force[0])

class RB_Dihed_F:
    def __init__(self, params, atomlist):
        self.params = params
        self.atomlist = atomlist

    def calc_forces(self, coords, force):
        for dihedral, params in self.params.items():
            atoms = dihedral.split('-')
            atomidx = [self.atomlist.index(atoms[0]), self.atomlist.index(atoms[1]), self.atomlist.index(atoms[2]), self.atomlist.index(atoms[3])]
            calc_rb_diheds(coords, atomidx, params, None, force[0])

class IMP_Dihed_F:
    def __init__(self, params, atomlist):
        self.params = params
        self.atomlist = atomlist

    def calc_forces(self, coords, force):
        for improper, [fconst, equ] in self.params.items():
            atoms = improper.split('-')
            atomidx = [self.atomlist.index(atoms[0]), self.atomlist.index(atoms[1]), self.atomlist.index(atoms[2]), self.atomlist.index(atoms[3])]
            calc_imp_diheds(coords, atomidx, equ, fconst, force[0])