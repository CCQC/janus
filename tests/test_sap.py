import pytest
from janus import qm_wrapper, mm_wrapper, qmmm
from copy import deepcopy
import numpy as np
import os

water = os.path.join(str('tests/files/test_openmm/water.pdb'))

psi4 = qm_wrapper.Psi4Wrapper()
openmm = mm_wrapper.OpenMMWrapper(sys_info=water,**{'md_ensemble':'NVT', 'return_info':[]})

openmm.initialize('Mechanical')
main_info_m = openmm.get_main_info()

sap_1 = qmmm.SAP(psi4, openmm, sys_info=water, aqmmm_param={'qmmm_param' : {'embedding_method' : 'Mechanical'}})
sap_2 = qmmm.SAP(psi4, openmm, sys_info=water, aqmmm_param={'qmmm_param' : {'embedding_method' : 'Mechanical'}})

sap_1.set_Rmin(2.6)
sap_1.set_Rmax(3.2)
sap_2.set_Rmin(2.6)
sap_2.set_Rmax(3.4)

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

    force1 = {0: np.array([  6.40880032e-05,   1.80584451e-02,   2.05000460e-02]), 
              3: np.array([ -5.69161234e-05,  -1.60375833e-02,  -1.82059526e-02]), 
              4: np.array([ -3.58593989e-06,  -1.01043090e-03,  -1.14704670e-03]), 
              5: np.array([ -3.58593989e-06,  -1.01043090e-03,  -1.14704670e-03])}
    force2 = {0: np.array([ 0.00144519,  0.42962668,  0.4845262 ]), 
              3: np.array([-0.00135455, -0.38167936, -0.43328451]),
              4: np.array([ -8.53419886e-05,  -2.40473028e-02,  -2.72986301e-02]), 
              5: np.array([ -8.53419886e-05,  -2.40473028e-02,  -2.72986301e-02]),
              6: np.array([  7.10846787e-05,   1.30801907e-04,   2.98005374e-03]),
              7: np.array([  4.47861467e-06,   8.24103519e-06,   1.87755120e-04]),
              8: np.array([  4.47861467e-06,   8.24103519e-06,   1.87755120e-04])}


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
    assert np.allclose(sap_1.systems[0]['qmmm_forces'][0], np.array([0.99983106, 0.952398,   0.94596195]))
    assert np.allclose(sap_2.systems[0]['qmmm_forces'][0], np.array([0.94494213,-15.49709206,-17.58769556]))
    assert np.allclose(sap_1.systems[0]['qmmm_forces'][1], np.ones((3))) 
    assert np.allclose(sap_2.systems[0]['qmmm_forces'][1], np.ones((3))) 
    assert len(sap_1.systems[0]['qmmm_forces']) == 6
    assert len(sap_2.systems[0]['qmmm_forces']) == 9

def test_run_qmmm():


    sap_1.run_qmmm(main_info_m, 'OpenMM')
    sap_2.run_qmmm(main_info_m, 'OpenMM')

    assert np.allclose(sap_1.systems[0]['qmmm_energy'],-0.007553840666010467)
    assert np.allclose(sap_2.systems[0]['qmmm_energy'],-0.007551817555662193)
    assert np.allclose(sap_1.systems[0]['qmmm_forces'][0], np.array([0.01119894, 0.04867413,-0.03788333]))
    assert np.allclose(sap_2.systems[0]['qmmm_forces'][0], np.array([0.01117866, 0.05045314,-0.03586097]))
    assert len(sap_1.systems[0]['qmmm_forces']) == 6
    assert len(sap_2.systems[0]['qmmm_forces']) == 9


