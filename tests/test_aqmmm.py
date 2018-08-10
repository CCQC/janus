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
    pass
#
#oxs.partition([0])
#for p in oxs.partitions:
#    print(p, oxs.partitions[p].qm_atoms)
#    
#oxs.save('qm', np.zeros((3,3)), 1.0)
#oxs.save('qm_bz', np.zeros((3,3)), 1.2)
#print(oxs.partitions['qm'].forces) 
#print(oxs.partitions['qm'].energy) 
#print(oxs.partitions['qm'].buffer_groups)
#
#
#s, d = oxs.compute_lamda_i(0.333580476259)
#print(s, d)
#
#s = oxs.get_switching_function(oxs.partitions['qm_bz'])
#a = oxs.get_info()
#                                                                                                   
