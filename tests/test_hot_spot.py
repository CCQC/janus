import pytest
from janus import qm_wrapper, mm_wrapper, qmmm
import numpy as np
import os

water = os.path.join(str('tests/files/test_openmm/water.pdb'))

psi4 = qm_wrapper.Psi4Wrapper()
openmm = mm_wrapper.OpenMMWrapper(sys_info=water,**{'md_ensemble':'NVT', 'return_info':[]})

openmm.initialize('Mechanical')
main_info_m = openmm.get_main_info()

hs =   qmmm.HotSpot(psi4, openmm, sys_info=water)
hs_0 = qmmm.HotSpot(psi4, openmm, sys_info=water)
hs_1 = qmmm.HotSpot(psi4, openmm, sys_info=water)
hs_2 = qmmm.HotSpot(psi4, openmm, sys_info=water)

hs_0.set_Rmin(2.6)
hs_0.set_Rmax(2.8)
hs_1.set_Rmin(2.6)
hs_1.set_Rmax(3.2)
hs_2.set_Rmin(2.6)
hs_2.set_Rmax(3.4)

def test_partition():

    hs.partition([0])
    hs_0.partition([0])
    hs_1.partition([0])
    hs_2.partition([0])

    assert np.allclose(hs.systems[0]['qm'].qm_atoms, np.array([0, 1, 2, 3, 4, 5, 6, 7, 8]))
    assert np.allclose(hs_0.systems[0]['qm'].qm_atoms, np.array([0, 1, 2]))
    assert np.allclose(hs_1.systems[0]['qm'].qm_atoms, np.array([0, 1, 2, 3, 4, 5]))
    assert np.allclose(hs_2.systems[0]['qm'].qm_atoms, np.array([0, 1, 2, 3, 4, 5, 6, 7, 8]))

def test_compute_lamda_i():

    lamda1 = hs.compute_lamda_i(3.0)
    lamda2 = hs.compute_lamda_i(3.8)
    lamda3 = hs.compute_lamda_i(4.0)
    lamda4 = hs.compute_lamda_i(4.5)
    lamda5 = hs.compute_lamda_i(5.0)

    assert lamda1[0] == 1
    assert lamda2[0] == 1
    assert np.allclose(lamda3[0], 0.8224337457798981)
    assert lamda4[0] == 0
    assert lamda5[0] == 0

def test_run_aqmmm():

    hs_0.systems[0]['qm'].qmmm_forces = {key: np.ones((1,3)) for key in range(3)}
    hs_1.systems[0]['qm'].qmmm_forces = {key: np.ones((1,3)) for key in range(6)}

    hs_0.run_aqmmm()
    hs_1.run_aqmmm()

    assert hs_0.systems[0]['qmmm_forces'] == hs_0.systems[0]['qm'].qmmm_forces
    assert np.allclose(hs_1.systems[0]['qmmm_forces'][0], np.ones((1,3))) 
    assert np.allclose(hs_1.systems[0]['qmmm_forces'][3], np.array([ 0.08217068,  0.08217068,  0.08217068]))
    
def test_run_qmmm():

    hs_0.run_qmmm(main_info_m, 'OpenMM')
    hs_1.run_qmmm(main_info_m, 'OpenMM')

    print(hs_0.systems[0]['qmmm_forces'][0])  
    print(hs_1.systems[0]['qmmm_forces'][0])
    print(hs_1.systems[0]['qmmm_forces'][3])

    assert np.allclose(hs_0.systems[0]['qmmm_forces'][0], np.array([ 0.01119897, 0.04866929,-0.03788886]))
    assert np.allclose(hs_1.systems[0]['qmmm_forces'][0], np.array([ 0.00984826, 0.04972007,-0.03577803]))
    assert np.allclose(hs_1.systems[0]['qmmm_forces'][3], np.array([-5.05189234e-04, 4.98408392e-05, 5.00504091e-03]))

    assert hs_0.run_ID == 1
    assert hs_1.run_ID == 1
