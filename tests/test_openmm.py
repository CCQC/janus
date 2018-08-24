"""
Testing for the openmm_wrapper module
"""
import pytest
from janus import openmm_wrapper 
import simtk.unit as OM_unit
import numpy as np
import os

#ala_water_pdb_file = os.path.join(str('tests/files/test_openmm/ala_water.pdb'))
water_pdb_file = os.path.join(str('tests/files/test_openmm/water.pdb'))
ala_pdb_file = os.path.join(str('tests/files/test_openmm/ala_ala_ala.pdb'))

config_water = {"mm_pdb_file" : water_pdb_file,
          "embedding_method" : "Electrostatic"}

config_ala = {"mm_pdb_file" : ala_pdb_file,
          "embedding_method" : "Electrostatic"}

wrapper = openmm_wrapper.OpenMM_wrapper(config_water)
wrapper_ala = openmm_wrapper.OpenMM_wrapper(config_ala)
#openmm_mech = openmm_wrapper.OpenMM_wrapper(sys_mech)
#openmm_elec = openmm_wrapper.OpenMM_wrapper(sys_elec)
#openmm_ala_link = openmm_wrapper.OpenMM_wrapper(sys_ala_link)
#openmm_ala_link_2 = openmm_wrapper.OpenMM_wrapper(sys_ala_link)


def test_create_new_residue_template():

    mod = wrapper_ala.create_modeller(keep_qm=False, qm_atoms=[0,1,2,3])
    wrapper_ala.create_new_residue_template(mod.topology)

    assert wrapper_ala.forcefield._templates['Modified_ALA'].name == 'Modified_ALA' 

def test_set_charge_zero():

    sys1 = wrapper.create_openmm_system(wrapper.pdb.topology)
    sys2 = wrapper.create_openmm_system(wrapper.pdb.topology)

    wrapper.set_charge_zero(sys1)
    wrapper.set_charge_zero(sys2, link_atoms = [0])

    assert sys1.getForce(3).getParticleParameters(0)[0]/OM_unit.elementary_charge == 0.0
    assert sys1.getForce(3).getParticleParameters(1)[0]/OM_unit.elementary_charge == 0.0
    assert sys2.getForce(3).getParticleParameters(0)[0]/OM_unit.elementary_charge == 0.0
    assert sys2.getForce(3).getParticleParameters(1)[0]/OM_unit.elementary_charge == 0.417

def test_set_LJ_zero():

    sys = wrapper.create_openmm_system(wrapper.pdb.topology)

    wrapper.set_LJ_zero(sys)

    assert sys.getForce(3).getParticleParameters(0)[0]/OM_unit.elementary_charge == -0.834
    assert sys.getForce(3).getParticleParameters(0)[1]/OM_unit.nanometer == 0.0
    assert sys.getForce(3).getParticleParameters(0)[2]/OM_unit.kilojoule_per_mole == 0.0

def test_create_openmm_system():
    sys_1 = wrapper.create_openmm_system(wrapper.pdb.topology)
    sys_2 = wrapper.create_openmm_system(wrapper.pdb.topology, include_coulomb='all', initialize=True)
    sys_3 = wrapper.create_openmm_system(wrapper.pdb.topology, include_coulomb='only')

    assert sys_1.getNumForces() == 5
    assert sys_2.getNumForces() == 6
    assert sys_3.getNumForces() == 2

#def test_compute_mm():
#    pass

#def test_initialize():
#    pass
#
#def test_take_step():
#    pass
#
#def test_get_main_info():
#    pass
#

def test_create_modeller():
    mod1 = wrapper_ala.create_modeller(qm_atoms=[0,1,2,3], keep_qm=True)
    mod2 = wrapper_ala.create_modeller(qm_atoms=[0,1,2,3], keep_qm=False)
    atom1 = mod1.topology.getNumAtoms()
    atom2 = mod2.topology.getNumAtoms()

    assert atom1 == 4
    assert atom2 == 29

