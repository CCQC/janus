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

das_1 = qmmm.DAS(psi4, openmm, sys_info=water, aqmmm_param={'qmmm_param' : {'embedding_method' : 'Mechanical'}})
das_2 = qmmm.DAS(psi4, openmm, sys_info=water, aqmmm_param={'qmmm_param' : {'embedding_method' : 'Mechanical'}})
das_1.set_Rmin(2.6)
das_1.set_Rmax(3.2)
das_2.set_Rmin(2.6)
das_2.set_Rmax(3.4)

def test_get_combos():

    b1 = [1]
    b2 = [1,2]
    das_1.define_buffer_zone([0])
    das_2.define_buffer_zone([0])
    das_1.qm_residues = [0]
    das_2.qm_residues = [0]
    combo1, sigma1 = das_1.get_combos(b1)
    combo2, sigma2 = das_2.get_combos(b2)

    assert(len(combo1)) == len(b1)
    assert(len(combo2)) == len(b2)
    assert np.allclose(sigma1[0], 0.068937282788963161)
    assert np.allclose(sigma2[0], 0.29333176789549764)
    assert np.allclose(sigma2[1], 0.018284135252690212)

def test_compute_lamda_i():
    l, d = das_1.compute_lamda_i(3.0)

    assert np.allclose(0.7407407407407413, l)
    assert d is None

def test_partition():

    das_1.partition([0])
    das_2.partition([0])

    print(len(das_1.systems[0])) 
    print(len(das_2.systems[0]))

    assert len(das_1.systems[0]) == 2
    assert len(das_2.systems[0]) == 3


def test_run_aqmmm():

    das_1.systems[0]['qm'].qmmm_forces = {key: np.ones((3)) for key in range(3)}
    das_1.systems[0][0].qmmm_forces = {key: np.ones((3)) for key in range(6)}
    das_2.systems[0]['qm'].qmmm_forces = {key: np.ones((3)) for key in range(3)}
    das_2.systems[0][0].qmmm_forces = {key: np.ones((3)) for key in range(6)}
    das_2.systems[0][1].qmmm_forces = {key: np.ones((3)) for key in range(9)}

    das_1.get_zero_energy() 
    das_2.get_zero_energy()

    das_1.run_aqmmm()
    das_2.run_aqmmm()

    assert np.allclose(das_1.systems[0]['qmmm_energy'],80.0300166424)
    assert np.allclose(das_2.systems[0]['qmmm_energy'],99.6972889591)
    assert np.allclose(das_1.systems[0]['qmmm_forces'][0], np.array([0.99861371, 0.99861371, 0.99861371]))
    assert np.allclose(das_2.systems[0]['qmmm_forces'][0], np.array([1.00000021, 1.00000021, 1.00000021]))
    assert len(das_1.systems[0]['qmmm_forces']) == 6
    assert len(das_2.systems[0]['qmmm_forces']) == 9

def test_run_qmmm():

    das_1.run_qmmm(main_info_m, 'OpenMM')
    das_2.run_qmmm(main_info_m, 'OpenMM')

    assert np.allclose(das_1.systems[0]['qmmm_energy'],-0.00753616966636)
    assert np.allclose(das_2.systems[0]['qmmm_energy'],-0.00752291994858)
    assert np.allclose(das_1.systems[0]['qmmm_forces'][0], np.array([ 0.01109033, 0.04867425,-0.03769082]))
    assert np.allclose(das_2.systems[0]['qmmm_forces'][0], np.array([ 0.01077123, 0.0489786 ,-0.03723941]))
    assert len(das_1.systems[0]['qmmm_forces']) == 6
    assert len(das_2.systems[0]['qmmm_forces']) == 9


