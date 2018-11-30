import pytest
import mdtraj as md
from janus import partition
import numpy as np
import os

water = os.path.join(str('tests/files/test_openmm/water.pdb'))
traj = md.load(water)

dis = partition.DistancePartition(traj, traj.topology, 3.8, 4.5)
dis_0 = partition.DistancePartition(traj, traj.topology, 3.8, 4.5)
dis_1 = partition.DistancePartition(traj, traj.topology, 3.8, 4.5)
dis_2 = partition.DistancePartition(traj, traj.topology, 3.8, 4.5)

def test_set_Rmin():

    dis_0.set_Rmin(2.6)
    dis_1.set_Rmin(2.6)
    dis_2.set_Rmin(2.6)
    assert dis.get_Rmin() == 3.8
    assert dis_0.get_Rmin() == 2.6
    assert dis_1.get_Rmin() == 2.6
    assert dis_2.get_Rmin() == 2.6

def test_set_Rmax():

    dis_0.set_Rmax(2.8)
    dis_1.set_Rmax(3.2)
    dis_2.set_Rmax(3.4)
    assert dis.get_Rmax() == 4.5
    assert dis_0.get_Rmax() == 2.8
    assert dis_1.get_Rmax() == 3.2
    assert dis_2.get_Rmax() == 3.4

def test_edit_atoms():

    atom1 = dis.edit_atoms(atoms=[0,1,2,3,4], res_idx=1, remove=True)
    atom2 = dis.edit_atoms(atoms=[0,1,2,3,4], res_idx=1, add=True)

    assert np.allclose(np.array([0,1,2]), np.array(atom1))
    assert np.allclose(np.array([0,1,2,3,4,5]), np.array(atom2))

    
def test_define_buffer_zone():
    
    dis.define_buffer_zone([0])
    dis_0.define_buffer_zone([0])
    dis_1.define_buffer_zone([0])
    dis_2.define_buffer_zone([0])
    
    assert (dis.buffer_atoms == [8] and not dis.buffer_groups)
    assert (not dis_0.buffer_atoms and not dis_0.buffer_groups)
    assert (dis_1.buffer_atoms == [3] and dis_1.buffer_groups[1].atoms == [3, 4, 5])
    assert np.allclose(dis_2.buffer_atoms, np.array([3, 5, 6]))
    assert (dis_2.buffer_groups[1].atoms == [3, 4, 5] and dis_2.buffer_groups[2].atoms == [6, 7, 8])


def test_get_residue_info():

    res = dis.get_residue_info(0)
    res1 = dis.get_residue_info(1)

    assert np.allclose(np.array([0,1,2]), np.array(res.atoms))
    assert np.allclose(0.0655606189723, res.r_i)
    assert np.allclose(np.array([3,4,5]), np.array(res1.atoms))
    assert np.allclose(3.10273031189, res1.r_i)

def test_compute_COM():

    xyz, weight, ratio = dis_1.compute_COM(atoms=[0,1,2])
    com = np.array([0.011130575, .354230624, .588089475])
    
    assert weight == {0: 15.999, 1: 1.008, 2: 1.008}
    assert ratio == {0: 0.8880932556203164, 1: 0.055953372189841796, 2: 0.055953372189841796}
    assert np.allclose(xyz, com)

