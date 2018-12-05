import pytest
from janus import qm_wrapper, mm_wrapper, qmmm
from copy import deepcopy
import numpy as np
import os

water = os.path.join(str('tests/files/test_openmm/water.pdb'))

psi4 = qm_wrapper.Psi4Wrapper()
openmm = mm_wrapper.OpenMMWrapper(sys_info=water, **{'md_ensemble':'NVT', 'return_info':[]})

openmm.initialize('Mechanical')
main_info_m = openmm.get_main_info()

pap_1 = qmmm.PAP(hl_wrapper=psi4, ll_wrapper=openmm, sys_info=water, qmmm_param={'embedding_method' : 'Mechanical'}, Rmin=2.6, Rmax=3.2)
pap_2 = qmmm.PAP(hl_wrapper=psi4, ll_wrapper=openmm, sys_info=water, qmmm_param={'embedding_method' : 'Mechanical'}, Rmin=2.6, Rmax=3.4)

def test_get_combos():

    buffers = [0,1,2,3,4]
    buffers2 = [0,1]
    buffer_distance = {0:1.2, 1:0.4, 2:0.8, 3:1.5, 4:0.2}
    buffer_distance2 = {0:3, 1:1}

    combo1 = pap_1.get_combos(buffers)
    combo2 = pap_1.get_combos(buffers2)

    assert(len(combo1)) == 2**len(buffers) - 1
    assert(len(combo2)) == 2**len(buffers2) - 1


def test_find_configurations():

    pap_1.find_buffer_zone()
    pap_2.find_buffer_zone()

    pap_1.find_configurations()
    pap_2.find_configurations()

    assert len(pap_1.systems[0]) == 2
    assert len(pap_2.systems[0]) == 4

def test_compute_sf_gradient():

    f1 =  pap_1.compute_sf_gradient()
    f2 = pap_2.compute_sf_gradient()

    force1 = {0: np.array([ 0.01971489,  5.55517746,  6.30626794]), 
              3: np.array([-0.01750866, -4.93351563, -5.60055402]), 
              4: np.array([-0.00110311, -0.31083091, -0.35285696]),
              5: np.array([-0.00110311, -0.31083091, -0.35285696])}
    force2 = {0: np.array([  0.39625842,   2.68001022,  18.5492068 ]), 
              3: np.array([-0.00618909, -1.74393472, -1.97972426]),
              4: np.array([-0.00038994, -0.10987475, -0.12473042]),
              5: np.array([-0.00038994, -0.10987475, -0.12473042]),
              6: np.array([ -0.34572534,  -0.63616428, -14.49370119]), 
              7: np.array([-0.02178206, -0.04008085, -0.91316025]),
              8: np.array([-0.02178206, -0.04008085, -0.91316025])}

    for i, f in f1.items():
        assert np.allclose(f, force1[i])
    for i, f in f2.items():
        assert np.allclose(f, force2[i])

def test_run_aqmmm():

    pap_1.systems[0]['qm'].qmmm_forces = {key: np.ones((3)) for key in range(3)}
    pap_1.systems[0][0].qmmm_forces = {key: np.ones((3)) for key in range(6)}
    pap_2.systems[0]['qm'].qmmm_forces = {key: np.ones((3)) for key in range(3)}
    pap_2.systems[0][0].qmmm_forces = {key: np.ones((3)) for key in range(6)}
    pap_2.systems[0][1].qmmm_forces = {key: np.ones((3)) for key in range(6)}
    pap_2.systems[0][2].qmmm_forces = {key: np.ones((3)) for key in range(9)}

    pap_1.get_zero_energy() 
    pap_2.get_zero_energy()

    pap_1.run_aqmmm()
    pap_2.run_aqmmm()

    assert np.allclose(pap_1.systems[0]['qmmm_energy'], 77.43370419740786) 
    assert np.allclose(pap_2.systems[0]['qmmm_energy'], 95.51933428487493)
    assert np.allclose(pap_1.systems[0]['qmmm_forces'][0], np.array([1.05036505, 15.19164921, 17.11043806]))
    assert np.allclose(pap_2.systems[0]['qmmm_forces'][0], np.array([1.17854883, 32.57782478, 39.51293667]))
    assert np.allclose(pap_1.systems[0]['qmmm_forces'][1], np.ones((3))) 
    assert np.allclose(pap_2.systems[0]['qmmm_forces'][1], np.ones((3))) 
    assert len(pap_1.systems[0]['qmmm_forces']) == 6
    assert len(pap_2.systems[0]['qmmm_forces']) == 9


def test_run_qmmm():

    pap_1.systems[0]['qm'].qmmm_forces = {key: np.ones((3)) for key in range(3)}
    pap_1.systems[0][0].qmmm_forces = {key: np.ones((3)) for key in range(6)}
    pap_2.systems[0]['qm'].qmmm_forces = {key: np.ones((3)) for key in range(3)}
    pap_2.systems[0][0].qmmm_forces = {key: np.ones((3)) for key in range(6)}
    pap_2.systems[0][1].qmmm_forces = {key: np.ones((3)) for key in range(6)}
    pap_2.systems[0][2].qmmm_forces = {key: np.ones((3)) for key in range(9)}

    pap_1.run_qmmm(main_info_m, 'OpenMM')
    pap_2.run_qmmm(main_info_m, 'OpenMM')

    assert np.allclose(pap_1.systems[0]['qmmm_energy'],  -0.007550404996134019)
    assert np.allclose(pap_2.systems[0]['qmmm_energy'], -0.007526130891361977)
    assert np.allclose(pap_1.systems[0]['qmmm_forces'][0], np.array([0.01115458, 0.04872366,-0.03779692] ))
    assert np.allclose(pap_2.systems[0]['qmmm_forces'][0], np.array([0.01083326, 0.04899138,-0.03727567] ))
    assert len(pap_1.systems[0]['qmmm_forces']) == 6
    assert len(pap_2.systems[0]['qmmm_forces']) == 9

