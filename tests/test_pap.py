import pytest
from janus import pap, psi4_wrapper, openmm_wrapper, initializer
from copy import deepcopy
import numpy as np
import os

water = os.path.join(str('tests/files/test_openmm/water.pdb'))

param = {"system" : {"mm_pdb_file": water},
         "qmmm" : {"embedding_scheme" : "Electrostatic"}}

config = initializer.Initializer(param, as_file=False)
psi4 = psi4_wrapper.Psi4_wrapper(config.qm_param)
openmm = openmm_wrapper.OpenMM_wrapper(config.mm_param)

pap_1 = pap.PAP(config.aqmmm_param, psi4, openmm)
pap_2 = pap.PAP(config.aqmmm_param, psi4, openmm)

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

def test_run_aqmmm():
    pass
