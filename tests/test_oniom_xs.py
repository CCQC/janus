import pytest
from janus import qm_wrapper, mm_wrapper, qmmm
import numpy as np
import os

water = os.path.join(str('tests/files/test_openmm/water.pdb'))

psi4 = qm_wrapper.Psi4Wrapper()
openmm = mm_wrapper.OpenMMWrapper(sys_info=water,**{'md_ensemble':'NVT', 'return_info':[]})

openmm.initialize('Mechanical')
main_info_m = openmm.get_main_info()

oxs =   qmmm.OniomXS(psi4, openmm, sys_info=water)
oxs_0 = qmmm.OniomXS(psi4, openmm, sys_info=water)
oxs_1 = qmmm.OniomXS(psi4, openmm, sys_info=water)
oxs_2 = qmmm.OniomXS(psi4, openmm, sys_info=water)

def test_set_Rmin():

    oxs_0.set_Rmin(2.6)
    oxs_1.set_Rmin(2.6)
    oxs_2.set_Rmin(2.6)
    assert oxs.get_Rmin() == 3.8
    assert oxs_0.get_Rmin() == 2.6
    assert oxs_1.get_Rmin() == 2.6
    assert oxs_2.get_Rmin() == 2.6

def test_set_Rmax():

    oxs_0.set_Rmax(2.8)
    oxs_1.set_Rmax(3.2)
    oxs_2.set_Rmax(3.4)
    assert oxs.get_Rmax() == 4.5
    assert oxs_0.get_Rmax() == 2.8
    assert oxs_1.get_Rmax() == 3.2
    assert oxs_2.get_Rmax() == 3.4

def test_edit_atoms():

    atom1 = oxs.edit_atoms(atoms=[0,1,2,3,4], res_idx=1, remove=True)
    atom2 = oxs.edit_atoms(atoms=[0,1,2,3,4], res_idx=1, add=True)

    assert np.allclose(np.array([0,1,2]), np.array(atom1))
    assert np.allclose(np.array([0,1,2,3,4,5]), np.array(atom2))

    
def test_define_buffer_zone():
    
    oxs.define_buffer_zone([0])
    oxs_0.define_buffer_zone([0])
    oxs_1.define_buffer_zone([0])
    oxs_2.define_buffer_zone([0])
    
    assert (oxs.buffer_atoms == [8] and not oxs.buffer_groups)
    assert (not oxs_0.buffer_atoms and not oxs_0.buffer_groups)
    assert (oxs_1.buffer_atoms == [3] and oxs_1.buffer_groups[1].atoms == [3, 4, 5])
    assert np.allclose(oxs_2.buffer_atoms, np.array([3, 5, 6]))
    assert (oxs_2.buffer_groups[1].atoms == [3, 4, 5] and oxs_2.buffer_groups[2].atoms == [6, 7, 8])


def test_get_residue_info():

    res = oxs.get_residue_info(0)
    res1 = oxs.get_residue_info(1)

    assert np.allclose(np.array([0,1,2]), np.array(res.atoms))
    assert np.allclose(0.0655606189723, res.r_i)
    assert np.allclose(np.array([3,4,5]), np.array(res1.atoms))
    assert np.allclose(3.10273031189, res1.r_i)

def test_partition():

    print('oxs')
    oxs.partition([0])
    print('oxs0')
    oxs_0.partition([0])
    print('oxs1')
    oxs_1.partition([0])
    print('oxs2')
    oxs_2.partition([0])

    assert np.allclose(  oxs.systems[0]['qm'].qm_atoms, np.array([0, 1, 2, 3, 4, 5, 6, 7,8]))
    assert np.allclose(oxs_0.systems[0]['qm'].qm_atoms, np.array([0, 1, 2]))
    assert np.allclose(oxs_1.systems[0]['qm'].qm_atoms, np.array([0, 1, 2]))
    assert np.allclose(oxs_1.systems[0]['qm_bz'].qm_atoms, np.array([0, 1, 2, 3, 4, 5]))
    assert np.allclose(oxs_2.systems[0]['qm'].qm_atoms, np.array([0, 1, 2]))
    assert np.allclose(oxs_2.systems[0]['qm_bz'].qm_atoms, np.array([0, 1, 2, 3, 4, 5, 6, 7, 8]))

def test_compute_lamda_i():

    s, d = oxs_1.compute_lamda_i(3.0)
    assert (np.allclose(s, 0.20987654320987748) and np.allclose(d, -.8230452674897127))

def test_compute_COM():

    xyz, weight, ratio = oxs_1.compute_COM(atoms=[0,1,2])
    com = np.array([0.011130575, .354230624, .588089475])
    
    assert weight == {0: 15.999, 1: 1.008, 2: 1.008}
    assert ratio == {0: 0.8880932556203164, 1: 0.055953372189841796, 2: 0.055953372189841796}
    assert np.allclose(xyz, com)


def test_get_switching_function():

    s_2 = oxs_2.get_switching_function()
    assert np.allclose(s_2, 0.13708445855549423)

def test_compute_zero_energy():

    oxs_0.compute_zero_energy()
    oxs_1.compute_zero_energy()
    oxs_2.compute_zero_energy()

    assert np.allclose(oxs_0.qm_zero_energies['HOH'], -74.96598998934344 ) 
    assert np.allclose(oxs_0.mm_zero_energies['HOH'], 3.5887724974514855e-08)
    assert np.allclose(oxs_1.qm_zero_energies['HOH'], -74.96598998934344 )
    assert np.allclose(oxs_1.mm_zero_energies['HOH'], 3.5887724974514855e-08)
    assert np.allclose(oxs_2.qm_zero_energies['HOH'], -74.96598998934344 )
    assert np.allclose(oxs_2.mm_zero_energies['HOH'], 3.5887724974514855e-08)
    
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
    
 
