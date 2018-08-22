"""
Testing for the openmm_wrapper module
"""
import pytest
from janus import openmm_wrapper 
from janus import system
import numpy as np
import os

ala_water_pdb_file = os.path.join(str('tests/files/test_openmm/ala_water.pdb'))
water_pdb_file = os.path.join(str('tests/files/test_openmm/water.pdb'))
ala_pdb_file = os.path.join(str('tests/files/test_openmm/ala_ala_ala.pdb'))

qmmm_elec = {"embedding_method" : "Electrostatic"}

sys_ala = system.System(qm={"qm_atoms": [0,1,2,3,4,5,6]}, mm={'mm_pdb_file' : ala_water_pdb_file})
sys_mech = system.System(qm={"qm_atoms": [0,1,2]}, mm={'mm_pdb_file' : water_pdb_file})
sys_elec = system.System(qmmm=qmmm_elec,qm={"qm_atoms": [0,1,2]}, mm={'mm_pdb_file' : water_pdb_file})
sys_ala_link = system.System(qm={"qm_atoms": [0, 1, 2, 3]}, mm={'mm_pdb_file' : ala_pdb_file})

openmm_ala = openmm_wrapper.OpenMM_wrapper(sys_ala)
openmm_mech = openmm_wrapper.OpenMM_wrapper(sys_mech)
openmm_elec = openmm_wrapper.OpenMM_wrapper(sys_elec)
openmm_ala_link = openmm_wrapper.OpenMM_wrapper(sys_ala_link)
openmm_ala_link_2 = openmm_wrapper.OpenMM_wrapper(sys_ala_link)


#@pytest.mark.datafiles('tests/examples/test_openmm/water.pdb')
#def test_get_state_info(datafiles):
#    """
#    Function to test that get_openmm_energy is getting energy correctly.
#    """
#    sys, pdb = create_system(datafiles, 'water.pdb')
#    sim = janus.openmm_wrapper.create_openmm_simulation(sys, 
#                                                        pdb.topology,
#                                                        pdb.positions)
#    state  = janus.openmm_wrapper.get_state_info(sim)
#    energy = state['potential'] + state['kinetic']
#    assert np.allclose(energy, -0.01056289236)
#
#
#
#def test_get_second_subsys():
#def test_get_entire_sys():
#def test_get_boundary():


def get_state_info():
    pass
def create_openmm_simulation():
    pass
def test_initialize():
    pass

def test_take_step():
    pass

def test_get_main_info():
    pass

def test_compute_mm():
    pass

def test_create_openmm_system():
    pass

def set_charge_zero():
    pass

    
def test_create_new_residue_template():
    mm = openmm_ala_link.create_modeller(keep_qm = False)
    openmm_ala_link.create_new_residue_template(mm.topology)
    
    assert  openmm_ala_link._ff._templates['Modified_ALA'].name == 'Modified_ALA'


