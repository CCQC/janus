import pytest
from janus import ap, psi4_wrapper, openmm_wrapper
from copy import deepcopy
import numpy as np
import os

water = os.path.join(str('tests/files/test_openmm/water.pdb'))

config_sap = {"mm_pdb_file" : water, "aqmmm_scheme" : "SAP", "embedding_scheme" : "Electrostatic"}
config_pap = {"mm_pdb_file" : water, "aqmmm_scheme" : "PAP", "embedding_scheme" : "Electrostatic"}
psi4 = psi4_wrapper.Psi4_wrapper(config_sap)
openmm = openmm_wrapper.OpenMM_wrapper(config_sap)

sap_1 = ap.AP(config_sap, psi4, openmm)
sap_2 = ap.AP(config_sap, psi4, openmm)
pap_1 = ap.AP(config_pap, psi4, openmm)
pap_2 = ap.AP(config_pap, psi4, openmm)

sap_1.set_Rmin(0.26)
sap_1.set_Rmax(0.32)
sap_2.set_Rmin(0.26)
sap_2.set_Rmax(0.34)

pap_1.set_Rmin(0.26)
pap_1.set_Rmax(0.32)
pap_2.set_Rmin(0.26)
pap_2.set_Rmax(0.34)

def test_get_combos():

    buffers = [0,1,2,3,4]
    buffers2 = [0,1]
    buffer_distance = {0:1.2, 1:0.4, 2:0.8, 3:1.5, 4:0.2}
    buffer_distance2 = {0:3, 1:1}

    combo1 = sap_1.get_combos(buffers, buffer_distance)
    order1 = np.array(deepcopy(sap_1.sap_order))
    combo2 = sap_1.get_combos(buffers2, buffer_distance2)
    order2 = np.array(deepcopy(sap_1.sap_order))
    combo3 = pap_1.get_combos(buffers, buffer_distance)
    combo4 = pap_1.get_combos(buffers2, buffer_distance2)

    assert(len(combo1)) == len(buffers)
    assert(len(combo2)) == len(buffers2)
    assert(len(combo3)) == 2**len(buffers) - 1
    assert(len(combo4)) == 2**len(buffers2) - 1
    assert np.allclose(order1, np.array([4,1,2,0,3]))
    assert np.allclose(order2, np.array([1,0]))


def test_partition():

    sap_1.partition([0])
    sap_2.partition([0])
    pap_1.partition([0])
    pap_2.partition([0])

    assert len(sap_1.systems[0]) == 2
    assert len(sap_2.systems[0]) == 3
    assert len(pap_1.systems[0]) == 2
    assert len(pap_2.systems[0]) == 4

def test_get_sap_switching_functions():

    func1 = sap_1.get_sap_switching_functions()
    func2 = sap_2.get_sap_switching_functions()
    
    assert func1 == {1: [3.5669066455610295e-05]}
    assert func2 == {1:[0.019398444948607135],2:[9.163106551223009e-08]}
    assert func2 == {1:[0.019398444948607135],2:[9.163106551223009e-08]}

def test_run_aqmmm():
    pass
