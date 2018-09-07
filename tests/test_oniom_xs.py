import pytest
from janus import oniom_xs, psi4_wrapper, openmm_wrapper
import numpy as np
import os

water = os.path.join(str('tests/files/test_openmm/water.pdb'))

config = {"mm_pdb_file" : water}
psi4 = psi4_wrapper.Psi4_wrapper(config)
openmm = openmm_wrapper.OpenMM_wrapper(config)
oxs =   oniom_xs.ONIOM_XS(config, psi4, openmm)
oxs_0 = oniom_xs.ONIOM_XS(config, psi4, openmm)
oxs_1 = oniom_xs.ONIOM_XS(config, psi4, openmm)
oxs_2 = oniom_xs.ONIOM_XS(config, psi4, openmm)

def test_set_Rmin():

    oxs_0.set_Rmin(0.26)
    oxs_1.set_Rmin(0.26)
    oxs_2.set_Rmin(0.26)
    assert oxs.get_Rmin() == 0.38
    assert oxs_0.get_Rmin() == 0.26
    assert oxs_1.get_Rmin() == 0.26
    assert oxs_2.get_Rmin() == 0.26

def test_set_Rmax():

    oxs_0.set_Rmax(0.28)
    oxs_1.set_Rmax(0.32)
    oxs_2.set_Rmax(0.34)
    assert oxs.get_Rmax() == 0.45
    assert oxs_0.get_Rmax() == 0.28
    assert oxs_1.get_Rmax() == 0.32
    assert oxs_2.get_Rmax() == 0.34

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


def test_partition():

    oxs.partition([0])
    oxs_0.partition([0])
    oxs_1.partition([0])
    oxs_2.partition([0])

    assert np.allclose(  oxs.systems[0]['qm'].qm_atoms, np.array([0, 1, 2, 3, 4, 5, 6, 7,8]))
    assert np.allclose(oxs_0.systems[0]['qm'].qm_atoms, np.array([0, 1, 2]))
    assert np.allclose(oxs_1.systems[0]['qm'].qm_atoms, np.array([0, 1, 2]))
    assert np.allclose(oxs_1.systems[0]['qm_bz'].qm_atoms, np.array([0, 1, 2, 3, 4, 5]))
    assert np.allclose(oxs_2.systems[0]['qm'].qm_atoms, np.array([0, 1, 2]))
    assert np.allclose(oxs_2.systems[0]['qm_bz'].qm_atoms, np.array([0, 1, 2, 3, 4, 5, 6, 7, 8]))


def test_compute_lamda_i():

    s, d = oxs_1.compute_lamda_i(0.30)
    assert (np.allclose(s, 0.20987654320987748) and np.allclose(d, -82.30452674897127))

def test_compute_COM():

    xyz, weight, ratio = oxs_1.compute_COM(atoms=[0,1,2])
    com = np.array([0.011130575, .354230624, .588089475])
    
    assert weight == {0: 15.999, 1: 1.008, 2: 1.008}
    assert ratio == {0: 0.8880932556203164, 1: 0.055953372189841796, 2: 0.055953372189841796}
    assert np.allclose(xyz, com)

def test_get_buffer_info():

    oxs_1.get_buffer_info()
    oxs_2.get_buffer_info()

    assert np.allclose(oxs_1.buffer_groups[1].s_i  , 0.0329177832296379)   
    assert np.allclose(oxs_2.buffer_groups[1].s_i  , 0.26960388601830587)
    assert np.allclose(oxs_2.buffer_groups[2].s_i  , 0.004565031092682581)
    assert np.allclose(oxs_1.buffer_groups[1].d_s_i, -29.7335089104)
    assert np.allclose(oxs_2.buffer_groups[1].d_s_i, -65.9020143551)
    assert np.allclose(oxs_2.buffer_groups[2].d_s_i, -6.12352511059) 


    assert np.allclose(oxs_1.buffer_distance[1], 0.31027303118865379)
    assert np.allclose(oxs_2.buffer_distance[1], 0.31027303118865379)
    assert np.allclose(oxs_2.buffer_distance[2], 0.3335804762589481 )

def test_get_switching_function():

    s_2 = oxs_2.get_switching_function()
    assert np.allclose(s_2, 0.13708445855549423)

def test_run_qmmm():
    pass
    

def test_run_aqmmm():
    pass
#
#    oxs.save('qm', np.zeros((9,3)), 1.0)
#    forces = oxs.get_info()
#    
#    forces_1 = oxs_1.get_info()
#    #forces_2 = oxs_2.get_info()
#
#    assert (np.allclose(forces, np.zeros((9,3))) and oxs.energy == 1.0)
#    assert (forces_1 is None and oxs_1.energy == 1.0065835566459276)
#    #assert (forces_2 is None and oxs_2.energy == 1.0274168917110988)
