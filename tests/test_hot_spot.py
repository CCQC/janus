import pytest
from janus import hot_spot, psi4_wrapper, openmm_wrapper
import numpy as np
import os

water = os.path.join(str('tests/files/test_openmm/water.pdb'))

config = {"mm_pdb_file" : water}
psi4 = psi4_wrapper.Psi4_wrapper(config)
openmm = openmm_wrapper.OpenMM_wrapper(config)

hs =   hot_spot.HotSpot(config, psi4, openmm)
hs_0 = hot_spot.HotSpot(config, psi4, openmm)
hs_1 = hot_spot.HotSpot(config, psi4, openmm)
hs_2 = hot_spot.HotSpot(config, psi4, openmm)

hs_0.set_Rmin(0.26)
hs_0.set_Rmax(0.28)
hs_1.set_Rmin(0.26)
hs_1.set_Rmax(0.32)
hs_2.set_Rmin(0.26)
hs_2.set_Rmax(0.34)

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

    lamda1 = hs.compute_lamda_i(0.30)
    lamda2 = hs.compute_lamda_i(0.38)
    lamda3 = hs.compute_lamda_i(0.40)
    lamda4 = hs.compute_lamda_i(0.45)
    lamda5 = hs.compute_lamda_i(0.50)

    assert lamda1[0] == 1
    assert lamda2[0] == 1
    assert lamda3[0] == 0.8224337457798975 
    assert lamda4[0] == 0
    assert lamda5[0] == 0

def test_run_qmmm():
    pass
    
def test_run_aqmmm():
    pass
