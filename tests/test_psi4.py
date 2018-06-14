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


def test_get_energy():
    """
    Function to test get_psi4_energy is getting energy correctly
    """

    psi4_mech.get_energy()
    psi4_elec.get_energy()

    assert np.allclose(psi4_mech._energy, -149.9288270081589)
    assert np.allclose(psi4_elec._energy, -149.93278909825366)

def test_get_gradient():
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


def test_get_qm():

    mech = psi4_mech.get_qm()
    elec = psi4_elec.get_qm()

    assert np.allclose(mech['energy'], psi4_mech._energy)
    assert np.allclose(mech['gradients'], psi4_mech._gradient)
    assert np.allclose(elec['energy'],   psi4_elec._energy)
    assert np.allclose(elec['gradients'], psi4_elec._gradient)

def test_get_scf_charges():

    psi4_mech.get_scf_charges()
    psi4_elec.get_scf_charges()

    charge_mech = np.array([-0.3742676, 0.18533499, 0.18958095, -0.37355466, 0.18900526, 0.18390106])
    charge_elec = np.array([-0.38374907, 0.19459701, 0.1898844, -0.37350301, 0.19019059, 0.18258009])

    assert np.allclose(psi4_mech._charges, charge_mech)
    assert np.allclose(psi4_elec._charges, charge_elec)
    
def test_get_energy_and_charges():

    psi4_mp2.get_energy_and_charges()

    charge_mp2 = np.array([-0.35524337, 0.17608227, 0.17979322, -0.35448027, 0.17916542, 0.17468273])

    assert np.allclose(psi4_mp2._charges, charge_mp2)
    assert np.allclose(psi4_mp2._energy,  -149.9997933859008)


def test_get_external_charge():

    # test link atom
    charge = psi4_elec.get_external_charges(link=True)
    crg = np.array([[-0.834, 0.115 ,  0.31300001 , 6.14799976],
                        [0.417, 0.12  ,  0.241      , 5.19299984],
                        [0.417, 0.66  , -0.41499998 , 6.44599974]])

    # test RC 
    charge_rc = psi4_rc.get_external_charges(RC=True)

    crg_rc = np.array([[0.0823, 0.81399999558925629, 0.86099997162818909, 1.4949999749660492],
                        [-0.1825, 2.0569999516010284, -0.7720000296831131, 1.2890000641345978], 
                        [0.0603, 3.1360000371932983, -0.75199998915195465, 1.0320000350475311],
                        [0.0603, 1.9900000095367432, -0.64099997282028198, 2.3950000107288361],
                        [0.0603, 1.656000018119812, -1.7820000648498535, 1.062999963760376],
                        [0.5973, 1.956000030040741, 1.5790000557899475, 0.03599999938160181],
                        [-0.5679, 1.2189999967813492, 2.5249999761581421, -0.20099999383091927],
                        [0.011233333333333333, 1.60150003, 0.97700004, 0.33600002],
                        [0.011233333333333333, 1.65199999, -0.19850001, 0.96250005],
                        [0.011233333333333333, 1.03050001, 0.61799999, 1.06550001]])

    # test RC 
    charge_rcd = psi4_rc.get_external_charges(RCD=True)
    crg_rcd = np.array([[0.07106666666666667, 0.81399999999999995, 0.86099996999999995, 1.4949999700000001],
                        [-0.19373333333333334, 2.0569999499999998, -0.77200002999999995, 1.28900006],
                        [0.0603, 3.1360000399999999, -0.75199998999999995, 1.03200004],
                        [0.0603, 1.9900000099999999, -0.64099996999999997, 2.39500001],
                        [0.0603, 1.65600002, -1.7820000600000001, 1.06299996],
                        [0.5860666666666667, 1.95600003, 1.57900006, 0.035999999999999997],
                        [-0.5679, 1.2190000000000001, 2.52499998, -0.20099998999999999],
                        [0.022466666666666666, 1.601500025, 0.97700003499999999, 0.33600002000000001],
                        [0.022466666666666666, 1.6519999849999998, -0.19850000999999998, 0.96250005000000005],
                        [0.022466666666666666, 1.0305000099999999, 0.61799998999999994, 1.0655000050000001]])

    assert np.allclose(np.array(charge), crg)
    assert np.allclose(np.array(charge_rc), crg_rc)
    assert np.allclose(np.array(charge_rcd), crg_rcd)




def test_get_redistributed_positions():

    bonds = psi4_rc._system.boundary_info['bonds_to_mm']
    mm = psi4_rc._system.boundary_info['mm_index']
    pos = psi4_rc._system.entire_sys['positions']

    new_pos = np.array([[ 1.60150003,  0.97700004,  0.33600002],
                       [ 1.65199999, -0.19850001,  0.96250005],
                       [ 1.03050001,  0.61799999,  1.06550001]])

    position = psi4_elec.get_redistributed_positions(pos, bonds, mm)
    
    assert np.allclose(position, new_pos)


    

