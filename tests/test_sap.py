import pytest
from janus import qm_wrapper, mm_wrapper, qmmm
from copy import deepcopy
import numpy as np
import os

water = os.path.join(str('tests/files/test_openmm/water.pdb'))

psi4 = qm_wrapper.Psi4Wrapper()
openmm = mm_wrapper.OpenMMWrapper(sys_info=water)

openmm.initialize('Mechanical')
main_info_m = openmm.get_main_info()

sap_1 = qmmm.SAP(psi4, openmm, sys_info=water, aqmmm_param={'qmmm_param' : {'embedding_method' : 'electrostatic'}})
sap_2 = qmmm.SAP(psi4, openmm, sys_info=water, aqmmm_param={'qmmm_param' : {'embedding_method' : 'electrostatic'}})

sap_1.set_Rmin(0.26)
sap_1.set_Rmax(0.32)
sap_2.set_Rmin(0.26)
sap_2.set_Rmax(0.34)

def test_get_combos():

    buffers = [0,1,2,3,4]
    buffers2 = [0,1]
    buffer_distance = {0:1.2, 1:0.4, 2:0.8, 3:1.5, 4:0.2}
    buffer_distance2 = {0:3, 1:1}

    combo1 = sap_1.get_combos(buffers, buffer_distance)
    order1 = np.array(deepcopy(sap_1.sap_order))
    combo2 = sap_1.get_combos(buffers2, buffer_distance2)
    order2 = np.array(deepcopy(sap_1.sap_order))

    assert(len(combo1)) == len(buffers)
    assert(len(combo2)) == len(buffers2)
    assert np.allclose(order1, np.array([4,1,2,0,3]))
    assert np.allclose(order2, np.array([1,0]))


def test_partition():

    sap_1.partition([0])
    sap_2.partition([0])

    assert len(sap_1.systems[0]) == 2
    assert len(sap_2.systems[0]) == 3

def test_get_switching_functions():

    sap_1.get_switching_functions()
    sap_2.get_switching_functions()
    
    assert np.allclose(sap_1.buffer_groups[1].chi_i, 29.37871636203129)
    assert np.allclose(sap_2.buffer_groups[1].chi_i, 2.7217256971620434)
    assert np.allclose(sap_2.buffer_groups[2].chi_i, 220.81236605454026)
    assert np.allclose(sap_1.buffer_groups[1].d_phi_i_scaler, -3.5224397927679847e-06)
    assert np.allclose(sap_2.buffer_groups[1].d_phi_i_scaler, -0.015636653418654025)
    assert np.allclose(sap_2.buffer_groups[2].d_phi_i_scaler, -1.2393051001903936e-09)
    assert np.allclose(sap_1.buffer_groups[1].phi_i, 3.5669066455610295e-05)
    assert np.allclose(sap_2.buffer_groups[1].phi_i, 0.019398444948607135)
    assert np.allclose(sap_2.buffer_groups[2].phi_i, 9.163106551223009e-08)

def test_compute_sf_gradients():

    f1 =  sap_1.compute_sf_gradient()
    f2 = sap_2.compute_sf_gradient()

    force1 =  {0: np.array([ 0.0064088 ,  1.80584451,  2.0500046 ]), 
                3: np.array([-0.00569161, -1.60375833, -1.82059526]), 
                4: np.array([-0.00035859, -0.10104309, -0.11470467]), 
                5: np.array([-0.00035859, -0.10104309, -0.11470467])}
    force2 = {0: np.array([  0.14451921,  42.96266838,  48.45262028]), 
                3: np.array([ -0.13545501, -38.16793622, -43.32845066]), 
                4: np.array([-0.0085342 , -2.40473028, -2.72986301]), 
                5: np.array([-0.0085342 , -2.40473028, -2.72986301]), 
                6: np.array([ 0.00710847,  0.01308019,  0.29800537]), 
                7: np.array([ 0.00044786,  0.0008241 ,  0.01877551]), 
                8: np.array([ 0.00044786,  0.0008241 ,  0.01877551])}

    for i, f in f1.items():
        assert np.allclose(f, force1[i])
    for i, f in f2.items():
        assert np.allclose(f, force2[i])

def test_run_aqmmm():

    sap_1.systems[0]['qm'].qmmm_forces = {key: np.ones((3)) for key in range(3)}
    sap_1.systems[0][0].qmmm_forces = {key: np.ones((3)) for key in range(6)}
    sap_2.systems[0]['qm'].qmmm_forces = {key: np.ones((3)) for key in range(3)}
    sap_2.systems[0][0].qmmm_forces = {key: np.ones((3)) for key in range(6)}
    sap_2.systems[0][1].qmmm_forces = {key: np.ones((3)) for key in range(9)}

    sap_1.get_zero_energy() 
    sap_2.get_zero_energy()

    sap_1.run_aqmmm()
    sap_2.run_aqmmm()

    assert np.allclose(sap_1.systems[0]['qmmm_energy'],74.96866395622227)
    assert np.allclose(sap_2.systems[0]['qmmm_energy'],76.42022722434459)
    assert np.allclose(sap_1.systems[0]['qmmm_forces'][0], np.array([ 0.98310642, -3.76020004, -4.40380521]))
    assert np.allclose(sap_2.systems[0]['qmmm_forces'][0], np.array([-4.50578712, -1648.70920742, -1857.76955754]))
    assert np.allclose(sap_1.systems[0]['qmmm_forces'][1], np.ones((3))) 
    assert np.allclose(sap_2.systems[0]['qmmm_forces'][1], np.ones((3))) 
    assert len(sap_1.systems[0]['qmmm_forces']) == 6
    assert len(sap_2.systems[0]['qmmm_forces']) == 9

def test_run_qmmm():

    sap_1.run_qmmm(main_info_m, 'OpenMM')
    sap_2.run_qmmm(main_info_m, 'OpenMM')

    print(sap_1.systems[0]['qmmm_energy'])
    print(sap_2.systems[0]['qmmm_energy'])
    print(sap_1.systems[0]['qmmm_forces'][0])
    print(sap_2.systems[0]['qmmm_forces'][0])
    assert sap_1.systems[0]['qmmm_energy'] == -0.02926774770488531 
    assert sap_2.systems[0]['qmmm_energy'] == -0.03066988296157528 
    assert np.allclose(sap_1.systems[0]['qmmm_forces'][0], np.array([1.91476321e-02, 5.77438451e+01, 5.49196637e+01]))
    assert np.allclose(sap_2.systems[0]['qmmm_forces'][0], np.array([ -0.17353223, 117.9511188, 107.91283519 ]))
    assert len(sap_1.systems[0]['qmmm_forces']) == 9
    assert len(sap_2.systems[0]['qmmm_forces']) == 9


