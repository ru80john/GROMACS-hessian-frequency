import sys
from Hessian.read_chk_file import QM
from Hessian.md_file_info import Atom, Bond, Angle, Dihedral, Improper
from Hessian.hessian import fit_hessian
from Hessian.thread import rotating_spinner
from Hessian.frequencies import calc_qm_vs_md_frequencies
from Hessian.term_force import VDW_Q_F, Bond_F, Angle_F, RB_Dihed_F, IMP_Dihed_F
import threading
import os

compound = 'BEN' # Molecule
mypath = './' # The path can find the directory of molecule

# The scaling value of vibrational frequecy, 0.957 is used for wb97XD/6-311++G(d,p) 
# Other scaling value can be found on https://cccbdb.nist.gov/vibscalejustx.asp
vib_scaling = 0.957 

include_nonbonded = False # True: Include nonbonded interaction in calculation
include_dihedral = False # True: Include diheral term in calculation


def run(compound, mypath, vib_scaling, include_dihedral, include_nonbonded):
    folder_path = f"{mypath}{compound}/FREQUENCY"

    if os.path.exists(folder_path) and os.path.isdir(folder_path):
        pass
    else:
        os.mkdir(folder_path)
        print('Folder "FREQUENCY" created!')

    qm = QM(compound, mypath, vib_scaling)
    
    atom = Atom(compound, mypath)
    atom_list = atom.atomlist
    bond_fconst = Bond(compound, mypath).bond_list_fconst
    angle_fconst = Angle(compound, mypath).angle_list_fconst
    imp_fconst = Improper(compound, mypath, include_dihedral).improper_list_fconst

    if include_dihedral and include_nonbonded:
        atom_fconst = atom.non_bonded_pairs
        dih_fconst = Dihedral(compound, mypath).dihedral_list_fconst
        md = [VDW_Q_F(atom_fconst), Bond_F(bond_fconst, atom_list), Angle_F(angle_fconst, atom_list), RB_Dihed_F(dih_fconst, atom_list), IMP_Dihed_F(imp_fconst, atom_list)]

    elif not include_dihedral and include_nonbonded:
        atom_fconst = atom.non_bonded_pairs
        md = [VDW_Q_F(atom_fconst), Bond_F(bond_fconst, atom_list), Angle_F(angle_fconst, atom_list), IMP_Dihed_F(imp_fconst, atom_list)]

    elif include_dihedral and not include_nonbonded:
        dih_fconst = Dihedral(compound, mypath).dihedral_list_fconst
        md = [Bond_F(bond_fconst, atom_list), Angle_F(angle_fconst, atom_list), RB_Dihed_F(dih_fconst, atom_list), IMP_Dihed_F(imp_fconst, atom_list)]
        
    else:
        md = [Bond_F(bond_fconst, atom_list), Angle_F(angle_fconst, atom_list), IMP_Dihed_F(imp_fconst, atom_list)]
    
    md_hessian_1d = fit_hessian(qm, md)
    calc_qm_vs_md_frequencies(mypath, compound, qm, md_hessian_1d)





if __name__ == "__main__":
    done_event = threading.Event()
    spinner_thread = threading.Thread(
        target=rotating_spinner, args=(done_event,))
    spinner_thread.start()
    try:
        print("Main code is running...")
        #################################################
        run(compound, mypath, vib_scaling, include_dihedral, include_nonbonded)
        #################################################  
    finally:
        done_event.set()
        spinner_thread.join()
    
    print("\nCalculation Done!")
    print(f'Results are saved in "{mypath}{compound}/FREQUENCY"')
    
