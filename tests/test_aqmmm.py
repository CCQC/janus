import pytest
from janus import oniom_xs
import numpy as np
import os

water = os.path.join(str('tests/files/test_openmm/water.pdb'))
# numbers denote how many buffer groups are in each
oxs = oniom_xs.ONIOM_XS('distance',   water)
oxs_0 = oniom_xs.ONIOM_XS('distance', water)
oxs_1 = oniom_xs.ONIOM_XS('distance', water)
oxs_2 = oniom_xs.ONIOM_XS('distance', water)

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
    assert (oxs_1.buffer_atoms == [3] and oxs_1.buffer_groups == {1: [3, 4, 5]})
    assert (np.allclose(oxs_2.buffer_atoms, np.array([3, 5, 6])) and oxs_2.buffer_groups == {1: [3, 4, 5], 2: [6, 7, 8]})

def test_partition():

    oxs.partition([0])
    oxs_0.partition([0])
    oxs_1.partition([0])
    oxs_2.partition([0])

    assert np.allclose(oxs.partitions['qm'].qm_atoms, np.array([0, 1, 2, 3, 4, 5, 6, 7,8]))
    assert np.allclose(oxs_0.partitions['qm'].qm_atoms, np.array([0, 1, 2]))
    assert np.allclose(oxs_1.partitions['qm'].qm_atoms, np.array([0, 1, 2]))
    assert np.allclose(oxs_1.partitions['qm_bz'].qm_atoms, np.array([0, 1, 2, 3, 4, 5]))
    assert np.allclose(oxs_2.partitions['qm'].qm_atoms, np.array([0, 1, 2]))
    assert np.allclose(oxs_2.partitions['qm_bz'].qm_atoms, np.array([0, 1, 2, 3, 4, 5, 6, 7, 8]))


def test_save():
    
    oxs_1.save('qm', np.zeros((9,3)), 1.0)
    oxs_1.save('qm_bz', np.ones((9,3)), 1.2)

    assert np.allclose(oxs_1.partitions['qm'].forces, np.zeros((9,3)))
    assert oxs_1.partitions['qm'].energy == 1.0
    assert np.allclose(oxs_1.partitions['qm_bz'].forces, np.ones((9,3)))
    assert oxs_1.partitions['qm_bz'].energy == 1.2

def test_compute_lamda_i():

    s, d = oxs_1.compute_lamda_i(0.333580476259)
    assert (np.allclose(s, 1.1588880014834282) and np.allclose(d, 2.3113804921728383))

def test_switching_function():

    oxs_2.save('qm', np.zeros((9,3)), 1.0)
    oxs_2.save('qm_bz', np.ones((9,3)), 1.2)

    s_1 = oxs_1.get_switching_function(oxs_1.partitions['qm_bz'])
    s_2 = oxs_2.get_switching_function(oxs_2.partitions['qm_bz'])

    assert np.allclose(s_1, 0.9670822167703623)
    assert np.allclose(s_2, 0.862915541444506)

def test_get_info():

    oxs.save('qm', np.zeros((9,3)), 1.0)
    forces = oxs.get_info()
    
    forces_1 = oxs_1.get_info()
    #forces_2 = oxs_2.get_info()

    assert (np.allclose(forces, np.zeros((9,3))) and oxs.energy == 1.0)
    assert (forces_1 is None and oxs_1.energy == 1.0065835566459276)
    #assert (forces_2 is None and oxs_2.energy == 1.0274168917110988)
    
