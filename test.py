# %%
from os import listdir
import numpy as np
from ase.io import read
from forces import calc_bonds, calc_angles, calc_rb_diheds, calc_imp_diheds, calc_pairs
import math

def contains_chars(input_string, check_chars):
    return all(chars in input_string for chars in check_chars)

class Rtp:
    def __init__(self, compound, root_path):
        self.compound = compound
        self.mypath = root_path

    def pdbinfo(self):
        return read(self.mypath + self.compound + '/STRUCTURE/' + self.compound + '.pdb')
    
    @property
    def rtplist(self):
        with open(self.mypath + self.compound + '/STRUCTURE/' + self.compound + '.pdb', 'r') as f:
            lines = f.readlines()
            rtp = []
            for i in lines:
                if 'HETATM' in i:
                    tmp = i.split()[3]
                    
                    # Determine the rtp file name
                    # rtp name in pdb file is different from rtp file name
                    ############################################
                    try:
                        with open('./' + self.compound + '/FORCEFIELD/' + self.compound + '.ff/' + tmp + '.rtp', 'r') as f:
                            pass
                    except FileNotFoundError:
                        mypath = './' + self.compound + '/FORCEFIELD/' + self.compound + '.ff/'
                        files = listdir(mypath)
                        temp = []
                        for i in files:
                            if contains_chars(i, tmp + '.rtp'):
                                temp.append(i)
                        if len(temp) >= 2:
                            tmp_tmp = []
                            for i in temp:
                                if len(i) >= 8:
                                    tmp_tmp.append(i)
                            if len(tmp_tmp) >= 2:
                                for i in tmp_tmp:
                                    if i[:2] == tmp[:2]:
                                        tmp = i[:-4]
                            else:
                                tmp = tmp_tmp[0][:-4]
                        else:
                            tmp = temp[0][:-4]
                    ############################################

                    if tmp not in rtp:
                        rtp.append(tmp)
        return rtp

class Atom(Rtp):
    fudgeLJ = 0.5
    fudgeQQ = 0.5
    comb_rule = 3

    @property
    def atomlist(self):
        atom_type = {}
        atom_list = []
        with open(self.mypath + self.compound + '/STRUCTURE/' + self.compound + '.pdb', 'r') as f:
            lines = f.readlines()
            for i in lines:
                if 'HETATM' in i:
                    tmp = i.split()
                    if atom_type.get(tmp[-1], 0) == 0:
                        atom_type[tmp[-1]] = 1
                        atom_list.append(f'{tmp[-1]}1')
                    else:
                        atom_type[tmp[-1]] += 1
                        atom_list.append(f'{tmp[-1]}{atom_type[tmp[-1]]}')
        return atom_list
    
    @property
    def coordinates(self):
        return self.pdbinfo().get_positions()
    
    @property
    def atomname(self):
        atom_name = {}
        idx = 0
        with open(self.mypath + self.compound + '/STRUCTURE/' + self.compound + '.pdb', 'r') as f:
            lines = f.readlines()
            for i in lines:
                if 'HETATM' in i:
                    tmp = i.split()
                    atom_name[self.atomlist[idx]] = tmp[2]
                    idx += 1
                elif 'END' in i:
                    break
        return atom_name
    @property
    def atomcharge(self):
        atom_charge = {}
        idx = 0
        for i in self.rtplist:
            with open(self.mypath + self.compound + '/FORCEFIELD/' + self.compound + '.ff/' + i + '.rtp', 'r') as f:
                lines = f.readlines()
                for index, i in enumerate(lines):
                    if 'atoms' in i:
                        start = index + 1
                        break
                for i in lines[start:]:
                    tmp = i.split()
                    try:
                        float(tmp[2])
                    except ValueError:
                        continue
                    except IndexError:
                        break
                    atom_charge[self.atomlist[idx]] = tmp[2]
                    idx += 1
        return atom_charge
    
    @property
    def atomtype(self):
        atom_type = {}
        idx = 0
        for i in self.rtplist:
            with open(self.mypath + self.compound + '/FORCEFIELD/' + self.compound + '.ff/' + i + '.rtp', 'r') as f:
                lines = f.readlines()
                for index, i in enumerate(lines):
                    if 'atoms' in i:
                        start = index + 1
                        break
                for i in lines[start:]:
                    tmp = i.split()
                    try:
                        float(tmp[2])
                    except ValueError:
                        continue
                    except IndexError:
                        break
                    atom_type[self.atomlist[idx]] = tmp[1]
                    idx += 1
        return atom_type
    
    @property
    def vdw_radius_depth(self):
        atomtype_radius_depth = {}
        atom_sig_eps = {}
        with open(self.mypath + self.compound + '/FORCEFIELD/' + self.compound + '.ff/' + 'forcefield.itp', 'r') as f:
            lines = f.readlines()
            for index, i in enumerate(lines):
                if 'atomtypes' in i:
                    atom_idx = index + 1
                    break
            for i in lines[atom_idx:]:
                tmp = i.split()
                try:
                    float(tmp[1])
                except ValueError:
                    continue
                except IndexError:
                    break
                atomtype_radius_depth[tmp[0]] = [float(tmp[4]), float(tmp[5])]

            for i in self.atomlist:
                atom_sig_eps[i] = atomtype_radius_depth[self.atomtype[i]]
             
        return atom_sig_eps
    @property
    def non_bonded_pairs(self):
        pair_sig_eps_q = {}
        with open(self.mypath + self.compound + '/FORCEFIELD/' + self.compound + '.top', 'r') as f:
            lines = f.readlines()
        
        for index, i in enumerate(lines):
            if 'pairs' in i:
                start = index + 2
            
            elif 'angles' in i:
                end = index - 1
                break
        for i in lines[start:end]:
            atom_pair = i.split()
            atom1 = self.atomlist[int(atom_pair[0]) - 1] 
            atom2 = self.atomlist[int(atom_pair[1]) - 1]

            comb_sig_eps = self.use_combination_rule(self.vdw_radius_depth[atom1], self.vdw_radius_depth[atom2], self.comb_rule)
            c6, c12 = self.sig_eps_to_c6_c12(comb_sig_eps, self.comb_rule)
            qq = self.charge_interact(float(self.atomcharge[atom1]), float(self.atomcharge[atom2]))
            
            if self.check_same_atp(atom1, atom2):
                sig_eps_q = np.array([c6, c12, qq])
                pair_sig_eps_q[[int(atom_pair[0]) - 1, int(atom_pair[1]) - 1]] = sig_eps_q
            else:
                c6 = c6 * self.fudgeLJ
                c12 = c12 * self.fudgeLJ
                qq = qq * self.fudgeQQ
                sig_eps_q = np.array([c6, c12, qq])
                pair_sig_eps_q[[int(atom_pair[0]) - 1, int(atom_pair[1]) - 1]] = sig_eps_q
        return pair_sig_eps_q  

    def calc_forces(self, force):
        for pair, params in self.non_bonded_pairs.items():
            force[0] = calc_pairs(self.coordinates, pair, params, force[0])
        return force 
    
    def check_same_atp(self, atom1, atom2):
        if self.atomtype[atom1] == self.atomtype[atom2]:
            return True
        else:
            return False

            
    @staticmethod
    def sig_eps_to_c6_c12(params, comb_rule):
        if comb_rule == 1:
            c6 = params[0] * 1e6
            c12 = params[1] * 1e12
        else:
            sigma = params[0] * 10
            epsilon = params[1]
            sigma6 = sigma**6
            c6 = 4 * epsilon * sigma6
            c12 = c6 * sigma6
        return c6, c12
    
    @staticmethod
    def use_combination_rule(param1, param2, comb_rule):
        b = (param1[1] * param2[1])**0.5
        if comb_rule in [1, 3]:
            a = (param1[0] * param2[0])**0.5
        else:
            a = (param1[0] + param2[0]) / 2
        return a, b
    @staticmethod
    def charge_interact(q1, q2):
        from ase.units import _eps0, kJ, mol, J, m
        inv_eps0 = 1/(4*np.pi*_eps0) * m / J / kJ * mol
        qq = q1 *q2 * inv_eps0
        return qq


    def gaussian_format(self):
        atom_atomtype = []
        file = read(self.mypath + self.compound + '/STRUCTURE/' + self.compound + '.pdb')
        atom_position = file.get_positions()
        atom_symbol = file.get_chemical_symbols()
        for idx, [x, y, z] in enumerate(atom_position):
            atom_atomtype.append(f'{atom_symbol[idx]}-{self.atomlist[idx]}  \t{x:>10}{y:>10}{z:>10}\n')
        
        atom_atomtype.append('\n')
        return atom_atomtype

    def eff_charge(self):
        atom_charge_format = []
        for atom, value in self.atomcharge.items():
            atom_charge_format.append(f'EffChg  {atom:<4}{value:>10}\n')
        return atom_charge_format

    def vdw(self):
        atom_vdw_format = []
        atom_vdw_format.append('NonBon 6 1, 0 0 0.0 0.0 1.0 0.0 0.0 1.0\n')
        for atom, [radius, depth] in self.vdw_radius_depth.items():
            atom_vdw_format.append(f'VDW {atom:<4}{radius:>10.6f}{depth:>10.6f}\n')
        return atom_vdw_format

class Bond(Atom):
    @property
    def bond_list(self):
        bondlist = []
        with open(self.mypath + self.compound + '/FORCEFIELD/' + self.compound + '.top', 'r') as f:
            lines = f.readlines()
            for index, i in enumerate(lines):
                if 'bonds' in i:
                    start = index + 2
                
                elif 'pairs' in i:
                    end = index - 1
                    break
            for i in lines[start:end]:
                tmp = i.split()
                bondlist.append('-'.join([self.atomlist[int(tmp[0]) - 1], self.atomlist[int(tmp[1]) - 1]]))
        return bondlist
    @property
    def bondtype_fconst(self):
        bondtypefconst = {}
        with open(self.mypath + self.compound + '/FORCEFIELD/' + self.compound + '.ff/' + 'forcefield.itp', 'r') as f:
            lines = f.readlines()
            for index, i in enumerate(lines):
                if 'bondtypes' in i:
                    start = index + 1
                elif 'angletypes' in i:
                    end = index
            for i in lines[start:end]:
                tmp = i.split()
                try:
                    float(tmp[4])
                except ValueError:
                    if ';' in tmp[4]:
                        cut = tmp[4].find(';')
                        tmp[4] = tmp[4][:cut]
                    else:
                        continue
                except IndexError:
                    continue
                if tmp[0][0] == ';':
                    continue
                bondtypefconst['-'.join([tmp[0], tmp[1]])] = [float(tmp[4])/100, float(tmp[3])*10]
                bondtypefconst['-'.join([tmp[1], tmp[0]])] = [float(tmp[4])/100, float(tmp[3])*10]
        return bondtypefconst
    
    @property
    def bond_list_fconst(self):
        bond_fconst = {}
        for i in self.bond_list:
            tmp = i.split('-')
            bond_fconst[i] = self.bondtype_fconst['-'.join([self.atomtype[tmp[0]], self.atomtype[tmp[1]]])]
        return bond_fconst
    
    def calc_forces(self, force):
        for bond, [fconst, b_length] in self.bond_list_fconst.items():
            atoms = bond.split('-')
            atomidx = [self.atomlist.index(atoms[0]), self.atomlist.index(atoms[1])]
            force[0] = calc_bonds(self.coordinates, atomidx, b_length, fconst, force[0])
        return force           
    
    def bond_fconst_format(self):
        bond_format = []
        for bond, [fconst, b_length] in self.bond_list_fconst.items():
            tmp = bond.split('-')
            bond_format.append(f'HrmStr1 {tmp[0]:<4}{tmp[1]:<4}{fconst:>6.1f}{b_length:>6.2f}\n')
        return bond_format

class Angle(Atom):
    '''
    HrmBnd1: Harmonic bend (Amber 1): ForceC*(θ - θeq)^2

    HrmBnd1 Atom-type1 Atom-type2 Atom-type3 ForceC θeq

    ForceC	: Force constant (in kcal/(mol*rad^2))
    θeq	: Equilibrium angle
    '''
    @property
    def angle_list(self):
        anglelist = []
        with open(self.mypath + self.compound + '/FORCEFIELD/' + self.compound + '.top', 'r') as f:
            lines = f.readlines()
            for index, i in enumerate(lines):
                if 'angles' in i:
                    start = index + 2
                elif 'dihedrals' in i:
                    end = index - 1
                    break
            for i in lines[start:end]:
                tmp = i.split()
                anglelist.append('-'.join([self.atomlist[int(tmp[0]) - 1], self.atomlist[int(tmp[1]) - 1], self.atomlist[int(tmp[2]) - 1]]))
        return anglelist
    
    @property
    def angletype_fconst(self):
        angletypefconst = {}
        with open(self.mypath + self.compound + '/FORCEFIELD/' + self.compound + '.ff/' + 'forcefield.itp', 'r') as f:
            lines = f.readlines()
            for index, i in enumerate(lines):
                if 'angletypes' in i:
                    start = index + 1
                elif 'dihedraltypes' in i:
                    end = index
            for i in lines[start:end]:
                tmp = i.split()
                try:
                    float(tmp[5])
                except ValueError:
                    if ';' in tmp[5]:
                        cut = tmp[5].find(';')
                        tmp[5] = tmp[5][:cut]
                    else:
                        continue
                except IndexError:
                    continue
                if tmp[0][0] == ';':
                    continue
                angletypefconst['-'.join([tmp[0], tmp[1], tmp[2]])] = [float(tmp[5]), math.radians(float(tmp[4]))]
                angletypefconst['-'.join([tmp[2], tmp[1], tmp[0]])] = [float(tmp[5]), math.radians(float(tmp[4]))]
        return angletypefconst

    @property
    def angle_list_fconst(self):
        angle_fconst = {}
        for i in self.angle_list:
            tmp = i.split('-')
            angle_fconst[i] = self.angletype_fconst['-'.join([self.atomtype[tmp[0]], self.atomtype[tmp[1]], self.atomtype[tmp[2]]])]
        return angle_fconst
    
    def calc_forces(self, force):
        for angle, [fconst, ang] in self.angle_list_fconst.items():
            atoms = angle.split('-')
            atomidx = [self.atomlist.index(atoms[0]), self.atomlist.index(atoms[1]), self.atomlist.index(atoms[2])]
            force[0] = calc_angles(self.coordinates, atomidx, ang, fconst, force[0])
            
        return force

    def angle_fconst_format(self):
        angle_format = []
        for angle, [fconst, ang] in self.angle_list_fconst.items():
            tmp = angle.split('-')
            angle_format.append(f'HrmBnd1 {tmp[0]:<4}{tmp[1]:<4}{tmp[2]:<4}{fconst:>7.1f}{ang:>7.1f}\n')
        return angle_format
            
class Dihedral(Atom):
    '''
    AmbTrs: Amber torsion (Amber 1): Σi=1,4 (Magi*[1+cos(i * θ - POi)])/NPaths

    AmbTrs Atom-type1 A-type2 A-type3 A-type4 PO1 PO2 PO3 PO4 Mag1 Mag2 Mag3 Mag4 NPaths

    PO1-PO4	: Phase offsets
    Mag1-Mag4	: V/2 magnitudes
    NPaths	: Number of paths. When zero or less, determined on-the-fly.
    '''

    @property
    def dihedral_list(self):
        dihedrallist = []
        with open(self.mypath + self.compound + '/FORCEFIELD/' + self.compound + '.top', 'r') as f:
            lines = f.readlines()
            count = 1
            for index, i in enumerate(lines):
                if 'dihedrals' in i and count == 1:
                    start = index + 2
                    count -= 1
                elif 'dihedrals' in i and count == 0:
                    end = index - 1
                    break
            for i in lines[start:end]:
                tmp = i.split()
                dihedrallist.append('-'.join([self.atomlist[int(tmp[0]) - 1], self.atomlist[int(tmp[1]) - 1], self.atomlist[int(tmp[2]) - 1], self.atomlist[int(tmp[3]) - 1]]))
        return dihedrallist
    
    @staticmethod
    def ryckaert_bellemans(constant, angle):
        import numpy as np
        val = 0
        for n, i in enumerate(constant):
            val += i*(np.cos(np.radians(angle-180)) ** n)
        return val / 4.184
    
    @staticmethod
    def amber_torsion(x, mag1, mag2, mag3, mag4, po1, po2, po3, po4):
        return mag1*(1 + np.cos(np.radians(x - po1))) + mag2*(1 + np.cos(np.radians(2*x - po2))) + mag3*(1 + np.cos(np.radians(3*x - po3))) + mag4*(1 + np.cos(np.radians(4*x - po4)))


    def dihedral_fitting(self, constant):
        #from scipy.optimize import curve_fit
        from scipy.optimize import differential_evolution
        x = np.linspace(-180 , 180, 1000)
        v = self.ryckaert_bellemans(constant, x)
        #params, covariance = curve_fit(self.amber_torsion, x, v)
        bounds = [(-100, 100), (-100, 100), (-100, 100), (-100, 100), (0, 720), (0, 720), (0, 720), (0, 720)]
        weights = 1
        result = differential_evolution(lambda params: np.sum(weights * (v - self.amber_torsion(x, *params))**2), bounds)
        params = result.x
        mag1, mag2, mag3, mag4, po1, po2, po3, p04 = params
        npath = 1
        return [po1, po2, po3, p04, mag1, mag2, mag3, mag4, npath]

    @property
    def dihedraltype_fconst(self):
        dihedraltypefconst = {}
        with open(self.mypath + self.compound + '/FORCEFIELD/' + self.compound + '.ff/' + 'forcefield.itp', 'r') as f:
            lines = f.readlines()
            for index, i in enumerate(lines):
                if 'dihedraltypes' in i:
                    start = index + 1

            for i in lines[start:]:
                tmp = i.split()
                try:
                    float(tmp[10])
                except ValueError:
                    if ';' in tmp[10]:
                        cut = tmp[10].find(';')
                        tmp[10] = tmp[10][:cut]
                    else:
                        continue
                except IndexError:
                    continue
                if tmp[0][0] == ';':
                    continue
                if tmp[4] != '3':
                    continue
                fconst = [float(i) for i in tmp[5:]]
                #fconst = self.dihedral_fitting(const)
                #file = self.pdbinfo()
                #dihe_angle = file.get_dihedral(int(tmp)) 
                #fconst = [0, 180 - (dihe_angle - 360), 0, 0, 0.0, 10.6775, 0.0, 0.0, 1.0]               
                dihedraltypefconst['-'.join([tmp[0], tmp[1], tmp[2], tmp[3]])] = fconst
                dihedraltypefconst['-'.join([tmp[3], tmp[2], tmp[1], tmp[0]])] = fconst
        return dihedraltypefconst
        
    @property
    def dihedral_list_fconst(self):
        dihedral_fconst = {}
        file = self.pdbinfo()
        for i in self.dihedral_list:
            tmp = i.split('-')
            dihedral_fconst[i] = self.dihedraltype_fconst['-'.join([self.atomtype[tmp[0]], self.atomtype[tmp[1]], self.atomtype[tmp[2]], self.atomtype[tmp[3]]])]
            #dihe_angle = file.get_dihedral(self.atomlist.index(tmp[0]), self.atomlist.index(tmp[1]), self.atomlist.index(tmp[2]), self.atomlist.index(tmp[3]))
            #dihedral_fconst[i] = [0, 180 + (int(dihe_angle) * 2), 0, 0, 0.0, 1.0, 0.0, 0.0, 1.0]
        return dihedral_fconst
    
    def calc_forces(self, force):
        for dihedral, params in self.dihedral_list_fconst.items():
            atoms = dihedral.split('-')
            atomidx = [self.atomlist.index(atoms[0]), self.atomlist.index(atoms[1]), self.atomlist.index(atoms[2]), self.atomlist.index(atoms[3])]
            force[0] = calc_rb_diheds(self.coordinates, atomidx, params, None, force[0])
        return force
    
    def dihedral_fconst_format(self):
        dihedral_format = []
        for dihedral, [po1, po2, po3, po4, mag1, mag2, mag3, mag4, npath] in self.dihedral_list_fconst.items():
            tmp = dihedral.split('-')
            dihedral_format.append(f'AmbTrs {tmp[0]:<4}{tmp[1]:<4}{tmp[2]:<4}{tmp[3]:<4}{int(po1):>6}{int(po2):>6}{int(po3):>6}{int(po4):>6}{mag1:>8.4f}{mag2:>9.4f}{mag3:>9.4f}{mag4:>9.4f}{npath:>5.1f}\n')
        return dihedral_format

class Improper(Bond):
    '''
    ImpTrs: Improper torsion (Amber 1): Mag*[1+cos(Period*(θ - PO))]

    ImpTrs Atom-type1 Atom-type2 Atom-type3 Atom-type4 Mag PO Period

    Mag	: V/2 Magnitude
    PO	: Phase offset
    Period	: Periodicity
    '''
    def __init__(self, compound, root_path, dihedral = True):
        super().__init__(compound, root_path)
        self.dihedral = dihedral

    def check_imp_format(self, a1, a2, a3, a4):
        atom1 = self.atomlist[a1 - 1]
        atom2 = self.atomlist[a2 - 1]
        atom3 = self.atomlist[a3 - 1]
        atom4 = self.atomlist[a4 - 1]
        chk_imp_type = []
        atom = atom1
        for i in [atom2, atom3, atom4]:
            if f'{atom}-{i}' in self.bond_list or f'{i}-{atom}' in self.bond_list:
                chk_imp_type.append(1)
            else:
                chk_imp_type.append(0)
            atom = i
        return chk_imp_type

    @property
    def improper_list(self):
        improperlist = []
        with open(self.mypath + self.compound + '/FORCEFIELD/' + self.compound + '.top', 'r') as f:
            lines = f.readlines()
            if self.dihedral:
                count = 1
            else:
                count = 0
            for index, i in enumerate(lines):
                if 'dihedrals' in i and count == 1:
                    count -= 1
                elif 'dihedrals' in i and count == 0:
                    start = index + 2
                elif 'Include Position restraint file' in i:
                    end = index - 1
                    break
            for i in lines[start:end]:
                tmp = i.split()
                chk_imp_type = self.check_imp_format(int(tmp[0]),int(tmp[1]),int(tmp[2]),int(tmp[3]))
                '''
                if chk_imp_type == [1, 0, 0]:
                    improperlist.append('-'.join([self.atomlist[int(tmp[0]) - 1], self.atomlist[int(tmp[1]) - 1], self.atomlist[int(tmp[2]) - 1], self.atomlist[int(tmp[3]) - 1], 'Imp']))
                elif chk_imp_type == [1, 1, 1]:
                    improperlist.append('-'.join([self.atomlist[int(tmp[0]) - 1], self.atomlist[int(tmp[1]) - 1], self.atomlist[int(tmp[2]) - 1], self.atomlist[int(tmp[3]) - 1], 'Amb']))
                elif chk_imp_type == [0, 0, 1]:
                    improperlist.append('-'.join([self.atomlist[int(tmp[3]) - 1], self.atomlist[int(tmp[2]) - 1], self.atomlist[int(tmp[1]) - 1], self.atomlist[int(tmp[0]) - 1], 'Imp']))
                '''
                if chk_imp_type == [1, 0, 0]:
                    improperlist.append('-'.join([self.atomlist[int(tmp[0]) - 1], self.atomlist[int(tmp[1]) - 1], self.atomlist[int(tmp[2]) - 1], self.atomlist[int(tmp[3]) - 1], 'Imp']))
                elif chk_imp_type == [1, 1, 1]:
                    improperlist.append('-'.join([self.atomlist[int(tmp[0]) - 1], self.atomlist[int(tmp[1]) - 1], self.atomlist[int(tmp[2]) - 1], self.atomlist[int(tmp[3]) - 1], 'Pro']))
                elif chk_imp_type == [0, 0, 1]:
                    improperlist.append('-'.join([self.atomlist[int(tmp[3]) - 1], self.atomlist[int(tmp[2]) - 1], self.atomlist[int(tmp[1]) - 1], self.atomlist[int(tmp[0]) - 1], 'Imp']))
        return improperlist

    @staticmethod
    def improper_fitting(const):
        import numpy as np
        from scipy.optimize import curve_fit
        def fit_func(x, mag, period):
            return mag * (1 + np.cos(np.radians(period * (x - 180))))
        def harmonic(x, const):
            return 1/2 * const/4.184 * np.radians((x - 0))**2
        x = np.linspace(-60, 60, 1000)
        v = harmonic(x, const)
        params, covariance = curve_fit(fit_func, x, v)
        mag , period= params
        po = 0
        return [mag, po, period]
  
    
    
    @property
    def impropertype_fconst(self):
        impropertypefconst = {}
        with open(self.mypath + self.compound + '/FORCEFIELD/' + self.compound + '.ff/' + 'forcefield.itp', 'r') as f:
            lines = f.readlines()
            for index, i in enumerate(lines):
                if 'dihedraltypes' in i:
                    start = index + 1

            for i in lines[start:]:
                tmp = i.split()
                try:
                    float(tmp[6])
                except ValueError:
                    if ';' in tmp[6]:
                        cut = tmp[6].find(';')
                        tmp[6] = tmp[6][:cut]
                    else:
                        continue
                except IndexError:
                    continue
                if tmp[0][0] == ';':
                    continue
                if tmp[4] != '2':
                    continue

                fconst = [float(tmp[6]), float(tmp[5])]
                impropertypefconst['-'.join([tmp[0], tmp[1], tmp[2], tmp[3]])] = fconst        
                

                '''
                if float(tmp[6]) == 500.000:
                    #impropertype_fconst['-'.join([tmp[0], tmp[1], tmp[2], tmp[3], 'Imp'])] = [127.5687, 0, 1]
                    #impropertype_fconst['-'.join([tmp[0], tmp[2], tmp[1], tmp[3], 'Imp'])] = [127.5687, 0, 1]
                    # Gaussian不能接受3位數
                    impropertypefconst['-'.join([tmp[0], tmp[1], tmp[2], tmp[3], 'Imp'])] = [99.00, 0, 1]
                    impropertypefconst['-'.join([tmp[0], tmp[2], tmp[1], tmp[3], 'Imp'])] = [99.00, 0, 1]
                elif float(tmp[6]) == 500.400:
                    #impropertype_fconst['-'.join([tmp[0], tmp[1], tmp[2], tmp[3], 'Imp'])] = [127.6708, 0, 1]
                    #impropertype_fconst['-'.join([tmp[0], tmp[2], tmp[1], tmp[3], 'Imp'])] = [127.6708, 0, 1]
                    # Gaussian不能接受3位數
                    impropertypefconst['-'.join([tmp[0], tmp[1], tmp[2], tmp[3], 'Imp'])] = [99.00, 0, 1]
                    impropertypefconst['-'.join([tmp[0], tmp[2], tmp[1], tmp[3], 'Imp'])] = [99.00, 0, 1]                  
                elif float(tmp[6]) == 167.400:
                    if tmp[3].strip() == 'H':
                        impropertypefconst['-'.join([tmp[0], tmp[1], tmp[2], tmp[3], 'Imp'])] = [42.71, 0, 1]
                        impropertypefconst['-'.join([tmp[0], tmp[2], tmp[1], tmp[3], 'Imp'])] = [42.71, 0, 1]
                    else:
                        impropertypefconst['-'.join([tmp[0], tmp[1], tmp[2], tmp[3], 'Amb'])] = [0, 180, 0, 0, 0.0, 10.6775, 0.0, 0.0, 1.0]
                        impropertypefconst['-'.join([tmp[3], tmp[2], tmp[1], tmp[0], 'Amb'])] = [0, 180, 0, 0, 0.0, 10.6775, 0.0, 0.0, 1.0]
                else:
                    fconst = self.improper_fitting(float(tmp[6]))
                    if fconst[0] > 99.00:
                        impropertypefconst['-'.join([tmp[0], tmp[1], tmp[2], tmp[3], 'Imp'])] = [99.00, 0, 1]
                        impropertypefconst['-'.join([tmp[0], tmp[2], tmp[1], tmp[3], 'Imp'])] = [99.00, 0, 1]
                    else:
                        impropertypefconst['-'.join([tmp[0], tmp[1], tmp[2], tmp[3], 'Imp'])] = fconst
                        impropertypefconst['-'.join([tmp[0], tmp[2], tmp[1], tmp[3], 'Imp'])] = fconst                        
                '''
        return impropertypefconst
    
    @property
    def improper_list_fconst(self):
        improper_fconst = {}
        for i in self.improper_list:
            tmp = i.split('-')
            if tmp[4] == 'Imp':
                try:
                    fconst = self.impropertype_fconst['-'.join([self.atomtype[tmp[0]], self.atomtype[tmp[1]], self.atomtype[tmp[2]], self.atomtype[tmp[3]]])]
                except KeyError:
                    fconst = self.impropertype_fconst['-'.join([self.atomtype[tmp[0]], self.atomtype[tmp[2]], self.atomtype[tmp[1]], self.atomtype[tmp[3]]])]
            elif tmp[4] == 'Pro':
                try:
                    fconst = self.impropertype_fconst['-'.join([self.atomtype[tmp[0]], self.atomtype[tmp[1]], self.atomtype[tmp[2]], self.atomtype[tmp[3]]])]
                except KeyError:
                    fconst = self.impropertype_fconst['-'.join([self.atomtype[tmp[3]], self.atomtype[tmp[2]], self.atomtype[tmp[1]], self.atomtype[tmp[0]]])]
            
            #improper_fconst[i] = self.impropertype_fconst['-'.join([self.atomtype[tmp[0]], self.atomtype[tmp[1]], self.atomtype[tmp[2]], self.atomtype[tmp[3]], tmp[4]])]
            improper_fconst[i[:-4]] = fconst     
        return improper_fconst
    
    def calc_forces(self, force):
        for improper, [fconst, equ] in self.improper_list_fconst.items():
            atoms = improper.split('-')
            atomidx = [self.atomlist.index(atoms[0]), self.atomlist.index(atoms[1]), self.atomlist.index(atoms[2]), self.atomlist.index(atoms[3])]
            force[0] = calc_imp_diheds(self.coordinates, atomidx, equ, fconst, force[0])
        return force

    def improper_fconst_format(self):
        improper_format = []
        Amb_format = []
        for improper, fconst in self.improper_list_fconst.items():
            tmp = improper.split('-')
            if len(fconst) <= 3:
                mag, po, period = fconst
                improper_format.append(f'ImpTrs {tmp[2]:<4}{tmp[1]:<4}{tmp[0]:<4}{tmp[3]:<4}{mag:>7.2f}{po:>7.2f}{period:>5.1f}\n')
            else:
                po1, po2, po3, po4, mag1, mag2, mag3, mag4, npath = fconst
                Amb_format.append(f'AmbTrs {tmp[0]:<4}{tmp[1]:<4}{tmp[2]:<4}{tmp[3]:<4}{po1:>4}{po2:>6}{po3:>4}{po4:>4}{mag1:>8.4f}{mag2:>9.4f}{mag3:>9.4f}{mag4:>9.4f}{npath:>5.1f}\n')
            
        return Amb_format + improper_format


if __name__ == '__main__':
    force = np.zeros((1, 65, 3))
    a = Bond('CNP', './')
    c = Angle('CNP', './')
    d = Improper('CNP', './')
    print(d.improper_list_fconst)
    

