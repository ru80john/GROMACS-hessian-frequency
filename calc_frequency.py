# %%
import sys
sys.path.append('C:/Users/ru80j_iujm6o5/work/Hessian_v3/Hessian')
from read_chk_file import QM
from md_file_info import Atom, Bond, Angle, Dihedral, Improper
from hessian import fit_hessian
from thread import rotating_spinner
from frequencies import calc_qm_vs_md_frequencies
from term_force import VDW_Q_F, Bond_F, Angle_F, RB_Dihed_F, IMP_Dihed_F
import threading
import os

compound = 'BEN'
mypath = './'
vib_scaling = 0.957
include_nonbonded = False
have_dihedral = False


def run(compound, mypath, vib_scaling, have_dihedral, include_nonbonded):
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
    imp_fconst = Improper(compound, mypath, have_dihedral).improper_list_fconst

    if have_dihedral and include_nonbonded:
        atom_fconst = atom.non_bonded_pairs
        dih_fconst = Dihedral(compound, mypath).dihedral_list_fconst
        md = [VDW_Q_F(atom_fconst), Bond_F(bond_fconst, atom_list), Angle_F(angle_fconst, atom_list), RB_Dihed_F(dih_fconst, atom_list), IMP_Dihed_F(imp_fconst, atom_list)]

    elif not have_dihedral and include_nonbonded:
        atom_fconst = atom.non_bonded_pairs
        md = [VDW_Q_F(atom_fconst), Bond_F(bond_fconst, atom_list), Angle_F(angle_fconst, atom_list), IMP_Dihed_F(imp_fconst, atom_list)]

    elif have_dihedral and not include_nonbonded:
        dih_fconst = Dihedral(compound, mypath).dihedral_list_fconst
        md = [Bond_F(bond_fconst, atom_list), Angle_F(angle_fconst, atom_list), RB_Dihed_F(dih_fconst, atom_list), IMP_Dihed_F(imp_fconst, atom_list)]
        
    else:
        md = [Bond_F(bond_fconst, atom_list), Angle_F(angle_fconst, atom_list), IMP_Dihed_F(imp_fconst, atom_list)]
    
    md_hessian_1d = fit_hessian(qm, md)
    calc_qm_vs_md_frequencies(mypath, compound, qm, md_hessian_1d)





if __name__ == "__main__":
    # 創建一個事件對象，用於通知旋轉符號的線程停止
    done_event = threading.Event()
    # 創建一個線程顯示旋轉符號
    spinner_thread = threading.Thread(
        target=rotating_spinner, args=(done_event,))
    # 啟動旋轉符號的線程
    spinner_thread.start()
    try:
        print("Main code is running...")
        #################################################
        run(compound, mypath, vib_scaling, have_dihedral, include_nonbonded)
        #################################################  
    finally:
        # 通知旋轉符號的線程停止
        done_event.set()
        # 等待旋轉符號的線程結束
        spinner_thread.join()
    
    print("\nCalculation Done!")
    print(f'Results are saved in "{mypath}{compound}/FREQUENCY"')
    
