import scipy.optimize as optimize
import numpy as np


def fit_hessian(qm, mol):
    hessian, full_md_hessian_1d = [], []
    non_fit = []
    qm_hessian = np.copy(qm.hessian)

    print("Calculating the MD hessian matrix elements...")
    full_md_hessian = calc_hessian(qm, mol)

    count = 0
    for i in range(qm.n_atoms*3):
        for j in range(i+1):
            hes = (full_md_hessian[i, j] + full_md_hessian[j, i]) / 2
            if all([h == 0 for h in hes]) or np.abs(qm_hessian[count]) < 0.0001:
                qm_hessian = np.delete(qm_hessian, count)
                full_md_hessian_1d.append(np.zeros(1))
            else:
                count += 1
                hessian.append(hes[:])
                full_md_hessian_1d.append(hes[:])
                non_fit.append(hes[-1])

    fit = 1
    
    full_md_hessian_1d = np.sum(full_md_hessian_1d * fit, axis=1)

    return full_md_hessian_1d


def calc_hessian(qm, mol):
    """
    Scope:
    -----
    Perform displacements to calculate the MD hessian numerically.
    """
    qm_coords = np.copy(qm.coords)
    full_hessian = np.zeros((3*qm.n_atoms, 3*qm.n_atoms, 1))
    for a in range(qm.n_atoms):
        for xyz in range(3):
            qm_coords[a][xyz] += 0.003
            f_plus = calc_force(qm_coords, mol, qm)
            qm_coords[a][xyz] -= 0.006
            f_minus = calc_force(qm_coords, mol, qm)
            qm_coords[a][xyz] += 0.003
            diff = - (f_plus - f_minus) / 0.006
            full_hessian[a*3+xyz] = diff.reshape(1, 3*qm.n_atoms).T
            #print(f'  ({a*3 + xyz + 1}/{int(qm.n_atoms*3)}) Done!')
    return full_hessian


def calc_force(coords, mol, qm):
    """
    Scope:
    ------
    For each displacement, calculate the forces from all terms.

    """

    force = np.zeros((1, qm.n_atoms, 3))

    for term in mol:
        term.calc_forces(coords, force)
    '''
    with mol.terms.add_ignore(['dihedral/flexible']):
        for term in mol.terms:
            term.do_fitting(coords, force)
    '''
    return force