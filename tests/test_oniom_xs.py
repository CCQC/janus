import pytest
from janus import qm_wrapper, mm_wrapper, qmmm
import numpy as np
import os

water = os.path.join(str('tests/files/test_openmm/water.pdb'))

psi4 = qm_wrapper.Psi4Wrapper()
openmm = mm_wrapper.OpenMMWrapper(sys_info=water,**{'md_ensemble':'NVT', 'return_info':[]})

openmm.initialize('Mechanical')
main_info_m = openmm.get_main_info()

oxs =   qmmm.OniomXS(psi4, openmm, sys_info=water, Rmin=3.8, Rmax=4.5)
oxs_0 = qmmm.OniomXS(psi4, openmm, sys_info=water, Rmin=2.6, Rmax=2.8)
oxs_1 = qmmm.OniomXS(psi4, openmm, sys_info=water, Rmin=2.6, Rmax=3.2)
oxs_2 = qmmm.OniomXS(psi4, openmm, sys_info=water, Rmin=2.6, Rmax=3.4)

def test_find_buffer_zone():

    oxs.find_buffer_zone()
    oxs_0.find_buffer_zone()
    oxs_1.find_buffer_zone()
    oxs_2.find_buffer_zone()
    
    assert not oxs.buffer_groups
    assert not oxs_0.buffer_groups
    assert oxs_1.buffer_groups[1].atoms == [3, 4, 5]
    assert (oxs_2.buffer_groups[1].atoms == [3, 4, 5] and oxs_2.buffer_groups[2].atoms == [6, 7, 8])

def test_find_configurations():

    oxs.find_configurations()
    oxs_0.find_configurations()
    oxs_1.find_configurations()
    oxs_2.find_configurations()

    assert np.allclose(  oxs.systems[0]['qm'].qm_atoms, np.array([0, 1, 2, 3, 4, 5, 6, 7,8]))
    assert np.allclose(oxs_0.systems[0]['qm'].qm_atoms, np.array([0, 1, 2]))
    assert np.allclose(oxs_1.systems[0]['qm'].qm_atoms, np.array([0, 1, 2]))
    assert np.allclose(oxs_1.systems[0]['qm_bz'].qm_atoms, np.array([0, 1, 2, 3, 4, 5]))
    assert np.allclose(oxs_2.systems[0]['qm'].qm_atoms, np.array([0, 1, 2]))
    assert np.allclose(oxs_2.systems[0]['qm_bz'].qm_atoms, np.array([0, 1, 2, 3, 4, 5, 6, 7, 8]))


def test_compute_lamda_i():

    s, d = oxs_1.compute_lamda_i(3.0)
    assert (np.allclose(s, 0.20987654320987748) and np.allclose(d, -.8230452674897127))

def test_get_switching_function():

    s_2 = oxs_2.get_switching_function()
    assert np.allclose(s_2, 0.13708445855549423)

def test_compute_zero_energy():

    oxs_0.compute_zero_energy()
    oxs_1.compute_zero_energy()
    oxs_2.compute_zero_energy()

    assert np.allclose(oxs_0.qm_zero_energies['HOH'], -74.96598998934344 ) 
    assert np.allclose(oxs_0.mm_zero_energies['HOH'], 2.294265421796016e-08)
    assert np.allclose(oxs_1.qm_zero_energies['HOH'], -74.96598998934344 )
    assert np.allclose(oxs_1.mm_zero_energies['HOH'], 2.294265421796016e-08)
    assert np.allclose(oxs_2.qm_zero_energies['HOH'], -74.96598998934344 )
    assert np.allclose(oxs_2.mm_zero_energies['HOH'], 2.294265421796016e-08)
    
def test_get_zero_energy():

    oxs_0.get_zero_energy()
    oxs_1.get_zero_energy()
    oxs_2.get_zero_energy()

    assert np.allclose(oxs_0.systems[0]['qm'].qmmm_energy, 74.96598998934344) 
    assert np.allclose(oxs_1.systems[0]['qm'].qmmm_energy, 74.96598998934344)
    assert np.allclose(oxs_2.systems[0]['qm'].qmmm_energy, 74.96598998934344)
    assert np.allclose(oxs_1.systems[0]['qm_bz'].qmmm_energy, 74.96598998934344 * 2) 
    assert np.allclose(oxs_2.systems[0]['qm_bz'].qmmm_energy, 74.96598998934344 * 3)

def test_run_aqmmm():

    oxs_0.systems[0]['qm'].qmmm_forces = {key: np.ones((1,3)) for key in range(3)}
    oxs_1.systems[0]['qm'].qmmm_forces = {key: np.ones((1,3)) for key in range(3)}
    oxs_1.systems[0]['qm_bz'].qmmm_forces = {key: np.ones((1,3)) for key in range(6)}

    oxs_0.run_aqmmm()
    oxs_1.run_aqmmm()

    assert oxs_0.systems[0]['qmmm_energy'] == oxs_0.systems[0]['qm'].qmmm_energy
    assert oxs_0.systems[0]['qmmm_forces'] == oxs_0.systems[0]['qm'].qmmm_forces
    assert np.allclose(oxs_1.systems[0]['qmmm_energy'],77.43370419740786)
    assert np.allclose(oxs_1.systems[0]['qmmm_forces'][0], np.ones((1,3))) is False

def test_run_qmmm():

    oxs_0.run_qmmm(main_info_m, 'OpenMM')
    oxs_1.run_qmmm(main_info_m, 'OpenMM')

    assert np.allclose(oxs_0.systems[0]['qmmm_energy'], -0.007553844392873543)
    assert np.allclose(oxs_0.systems[0]['qmmm_forces'][0], np.array([ 0.01119897, 0.04866929,-0.03788886]))
    assert np.allclose(oxs_1.systems[0]['qmmm_energy'],-0.007550404996134019)
    assert np.allclose(oxs_1.systems[0]['qmmm_forces'][0], np.array([ 0.01115458,  0.04872366, -0.03779692]))
    assert oxs_0.run_ID == 1
    assert oxs_1.run_ID == 1
    
 
