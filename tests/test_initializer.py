import pytest
from janus import initializer
import mdtraj as md
import numpy as np
import os


water = os.path.join(str('tests/files/test_openmm/water.pdb'))
param1 = {"system" : {"mm_pdb_file" : water}}
param2 = {
            "system" : {"mm_pdb_file" : water},

            "aqmmm" : {"aqmmm_scheme" : "Hot-Spot",
                        "aqmmm_partition_scheme" : "distance",
                        "qm_center" : [0]},

            "qmmm" : {"embedding_method" : "Electrostatic",
                    "qm_program" : "Psi4",
                    "mm_program" : "OpenMM"}}

param3 = os.path.join(str('tests/files/test_initializer/input.json'))

init1 = initializer.Initializer(param1, as_file=False)
init2 = initializer.Initializer(param2, as_file=False)
init3 = initializer.Initializer(param3)

def test_initialize_wrappers():

    mm1, qmmm1 = init1.initialize_wrappers()
    mm2, qmmm2 = init2.initialize_wrappers()
    mm3, qmmm3 = init3.initialize_wrappers()

    assert qmmm1.class_type == 'ONIOM-XS'
    assert qmmm2.class_type == 'Hot-Spot'
    assert qmmm3.class_type == 'QMMM'
