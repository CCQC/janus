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

pap_1 = qmmm.PAP(psi4, openmm, sys_info=water, aqmmm_param={'qmmm_param' : {'embedding_method' : 'electrostatic'}})
pap_2 = qmmm.PAP(psi4, openmm, sys_info=water, aqmmm_param={'qmmm_param' : {'embedding_method' : 'electrostatic'}})

pap_1.set_Rmin(0.26)
pap_1.set_Rmax(0.32)
pap_2.set_Rmin(0.26)
pap_2.set_Rmax(0.34)

def test_get_combos():

    buffers = [0,1,2,3,4]
    buffers2 = [0,1]
    buffer_distance = {0:1.2, 1:0.4, 2:0.8, 3:1.5, 4:0.2}
    buffer_distance2 = {0:3, 1:1}

    combo1 = pap_1.get_combos(buffers)
    combo2 = pap_1.get_combos(buffers2)

    assert(len(combo1)) == 2**len(buffers) - 1
    assert(len(combo2)) == 2**len(buffers2) - 1


def test_partition():

    pap_1.partition([0])
    pap_2.partition([0])

    assert len(pap_1.systems[0]) == 2
    assert len(pap_2.systems[0]) == 4

def test_compute_sf_gradient():

    f1 =  pap_1.compute_sf_gradient()
    f2 = pap_2.compute_sf_gradient()

    force1 = {0: np.array([ 1.97148884,  555.51774573,  630.6267936 ]), 
              3: np.array([ -1.75086594, -493.35156336, -560.05540221]), 
              4: np.array([ -0.11031145, -31.08309118, -35.28569569]),
              5: np.array([ -0.11031145, -31.08309118, -35.28569569])}

    force2 =  {0:  np.array([   39.62584163,   268.0010218 ,  1854.92067971]),
               3: np.array([  -0.61890873, -174.39347167, -197.9724261 ]),
               4: np.array([ -0.03899369, -10.98747543, -12.47304241]),
               5: np.array([ -0.03899369, -10.98747543, -12.47304241]),
               6: np.array([  -34.57253397,   -63.61642829, -1449.37011926]),
               7: np.array([ -2.17820578,  -4.00808549, -91.31602477]),
               8: np.array([ -2.17820578,  -4.00808549, -91.31602477])}

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
    assert np.allclose(pap_1.systems[0]['qmmm_forces'][0], np.array([6.03650482,  1420.16491984,  1612.04380529]))
    assert np.allclose(pap_2.systems[0]['qmmm_forces'][0], np.array([   18.85488287,  3158.7824764 ,  3852.29366477]))
    assert np.allclose(pap_1.systems[0]['qmmm_forces'][1], np.ones((3))) 
    assert np.allclose(pap_2.systems[0]['qmmm_forces'][1], np.ones((3))) 
    assert len(pap_1.systems[0]['qmmm_forces']) == 6
    assert len(pap_2.systems[0]['qmmm_forces']) == 9


def test_run_qmmm():

    pap_1.run_qmmm(main_info_m, 'OpenMM')
    pap_2.run_qmmm(main_info_m, 'OpenMM')

    assert pap_1.systems[0]['qmmm_energy'] == -0.03977812308047926
    assert pap_2.systems[0]['qmmm_energy'] == -0.041810966263750866
    assert np.allclose(pap_1.systems[0]['qmmm_forces'][0], np.array([-0.03073274,-0.47845509,-1.07936284]))
    assert np.allclose(pap_2.systems[0]['qmmm_forces'][0], np.array([-0.01825466,-0.234095  ,-0.59877961 ]))
    assert len(pap_1.systems[0]['qmmm_forces']) == 9
    assert len(pap_2.systems[0]['qmmm_forces']) == 9
