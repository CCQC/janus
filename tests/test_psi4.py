"""
Testing for psi4_wrapper.py module
"""
import pytest
from janus import psi4_wrapper
from janus import system
import numpy as np
from copy import deepcopy

qmmm_elec = {"embedding_method" : "Electrostatic"}
QM = {"qm_atoms" : [0,1,2,3,4,5]}
qm_rc = {"qm_atoms" : [0,1,2,3]}

sys_mech = system.System(qm=QM)
sys_elec = system.System(qmmm=qmmm_elec, qm=QM)
sys_rc = system.System(qmmm=qmmm_elec, qm=qm_rc)

qm_mol = """O     0.123   3.593   5.841 
 H    -0.022   2.679   5.599 
 H     0.059   3.601   6.796 
 O     0.017   6.369   7.293 
 H    -0.561   5.928   6.669 
 H     0.695   6.771   6.749 
 """

sys_mech.qm_positions = qm_mol
sys_elec.qm_positions = qm_mol

sys_elec.second_subsys = {"charges" : [-0.834, 0.417, 0.417], "positions" : np.asarray([[ 0.115 ,  0.31300001 , 6.14799976],
                                                                                       [ 0.12  ,  0.241      , 5.19299984],
                                                                                       [ 0.66  , -0.41499998 , 6.44599974]])}


sys_rc.boundary_info['mm_index'] = 4
sys_rc.boundary_info['bonds_to_mm'] = [10, 6, 5]

sys_rc.entire_sys = {"charges" : [-0.4157, 0.2719, 0.1997, 0.1997, 0.0337, 0.0823, -0.1825, 0.0603, 0.0603, 0.0603, 0.5973, -0.5679],
                     "positions" : np.array([[ 0.024    ,  -0.103   ,  -0.101     ], 
                                            [ 0.027     , -1.13200001, -0.239     ],
                                            [-0.80499999,  0.163     ,  0.471     ],
                                            [-0.059     ,  0.38400002, -1.01899996],
                                            [ 1.24700002,  0.37500001,  0.63600004],
                                            [ 0.814     ,  0.86099997,  1.49499997],
                                            [ 2.05699995, -0.77200003,  1.28900006],
                                            [ 3.13600004, -0.75199999,  1.03200004],
                                            [ 1.99000001, -0.64099997,  2.39500001],
                                            [ 1.65600002, -1.78200006,  1.06299996],
                                            [ 1.95600003,  1.57900006,  0.036     ],
                                            [ 1.219     ,  2.52499998, -0.20099999]])}
          
          
psi4_mech = psi4_wrapper.Psi4_wrapper(sys_mech)
psi4_elec = psi4_wrapper.Psi4_wrapper(sys_elec)
psi4_rc = psi4_wrapper.Psi4_wrapper(sys_rc)
psi4_mp2 = deepcopy(psi4_mech) 
psi4_mp2._system.qm_method = 'mp2'


def test_compute_energy():
    """
    Function to test get_psi4_energy is getting energy correctly
    """

    psi4_mech.get_energy()
    psi4_elec.get_energy()

    assert np.allclose(psi4_mech._energy, -149.9288270081589)
    assert np.allclose(psi4_elec._energy, -149.93278909825366)

def test_compute_gradient():
    """
    Function to test get_psi4_gradient is getting the gradient correctly
    """

    psi4_mech.get_gradient()
    psi4_elec.get_gradient()

    gradient_mech = np.array([[-0.00995519, -0.04924799,  0.03606667],
                              [ 0.00685005,  0.0390004 ,  0.00096903],
                              [ 0.00308006,  0.00941025, -0.03776953],
                              [ 0.00633522, -0.00137517, -0.0604841 ],
                              [ 0.01545583,  0.01344479,  0.03199705],
                              [-0.02176597, -0.01123228,  0.02922087]])

    gradient_elec = np.array([[-0.00984083, -0.04908026,  0.03652675],
                              [ 0.00674015,  0.03945264, -0.00167509],
                              [ 0.00310424,  0.01065926, -0.036673  ],
                              [ 0.00656583, -0.00182659, -0.06047419],
                              [ 0.01528687,  0.01378971,  0.03179818],
                              [-0.02186328, -0.01109774,  0.02931768]])
  
    assert np.allclose(psi4_mech._gradient, gradient_mech)
    assert np.allclose(psi4_elec._gradient, gradient_elec)


def test_run_qm():

    mech = psi4_mech.get_qm(qm_positions=None)
    elec = psi4_elec.get_qm(qm_positions=None)

    assert np.allclose(mech['energy'], psi4_mech._energy)
    assert np.allclose(mech['gradients'], psi4_mech._gradient)
    assert np.allclose(elec['energy'],   psi4_elec._energy)
    assert np.allclose(elec['gradients'], psi4_elec._gradient)

def test_compute_scf_charges():

    psi4_mech.get_scf_charges()
    psi4_elec.get_scf_charges()

    charge_mech = np.array([-0.3742676, 0.18533499, 0.18958095, -0.37355466, 0.18900526, 0.18390106])
    charge_elec = np.array([-0.38374907, 0.19459701, 0.1898844, -0.37350301, 0.19019059, 0.18258009])

    assert np.allclose(psi4_mech._charges, charge_mech)
    assert np.allclose(psi4_elec._charges, charge_elec)
    
def test_compute_energy_and_charges():

    psi4_mp2.get_energy_and_charges()

    charge_mp2 = np.array([-0.35524337, 0.17608227, 0.17979322, -0.35448027, 0.17916542, 0.17468273])

    assert np.allclose(psi4_mp2._charges, charge_mp2)
    assert np.allclose(psi4_mp2._energy,  -149.9997933859008)


def test_build_qm_param():
    
    assert sys1.qm_param['basis'] == 'STO-3G'
    assert sys1.qm_param['scf_type'] == 'df' 
    assert sys1.qm_param['guess'] == 'sad' 
    assert sys1.qm_param['reference'] == 'rhf'
    assert sys1.qm_param['e_convergence'] == 1e-8
    assert sys1.qm_param['d_convergence'] == 1e-8 

    assert sys2.qm_param['basis'] == '6-31G'
    assert sys2.qm_param['scf_type'] == 'pk' 
    assert sys2.qm_param['guess'] == 'sad' 
    assert sys2.qm_param['reference'] == 'rhf'
    assert sys2.qm_param['e_convergence'] == 1e-9
    assert sys2.qm_param['d_convergence'] == 1e-8 
