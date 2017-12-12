"""
Testing for the openmm_wrapper module
"""
import pytest
import janus
import numpy as np
import os


def create_system(datafiles, filename):
    """
    Function for creating a Janus system with
    mm options used for testing.
    Tests on input_protein.pdb found in tests/examples/test_openmm
    """
    path = str(datafiles)
    mm_pdb_file = os.path.join(path, filename)
    pdb = janus.openmm_wrapper.create_openmm_pdb(mm_pdb_file)

    openmm_sys = janus.openmm_wrapper.create_openmm_system(pdb.topology)

    return openmm_sys, pdb


@pytest.mark.datafiles('tests/examples/test_openmm/water.pdb')
def test_get_openmm_energy(datafiles):
    """
    Function to test that get_openmm_energy is getting energy correctly.
    """
    sys, pdb = create_system(datafiles, 'water.pdb')
    sim = janus.openmm_wrapper.create_openmm_simulation(sys, 
                                                        pdb.topology,
                                                        pdb.positions)
    state  = janus.openmm_wrapper.get_state_info(sim)
    energy = state['potential'] + state['kinetic']
    assert np.allclose(energy, -27.732873913249932)



@pytest.mark.datafiles('tests/examples/test_openmm/ala_water.pdb')
def test_keep_residues(datafiles):
    """
    Function to test keep_residues function.
    """
    sys, pdb = create_system(datafiles, 'ala_water.pdb')
    mod_str = janus.openmm_wrapper.create_openmm_modeller(pdb)
    mod_int = janus.openmm_wrapper.create_openmm_modeller(pdb)
    janus.openmm_wrapper.keep_residues(mod_str, ['ALA'])
    janus.openmm_wrapper.keep_residues(mod_int, [0,1,2])
    res_str = mod_str.topology.getNumResidues()
    res_int = mod_int.topology.getNumResidues()
    assert res_str == 3 
    assert res_int == 3 

@pytest.mark.datafiles('tests/examples/test_openmm/ala_water.pdb')
def test_keep_atoms(datafiles):
    """
    Function to test keep_atoms function.
    """
    qm_atm = [0, 1, 2, 3, 4, 5, 6]
    sys, pdb = create_system(datafiles, 'ala_water.pdb')
    mod_str = janus.openmm_wrapper.create_openmm_modeller(pdb)
    mod_int = janus.openmm_wrapper.create_openmm_modeller(pdb)
    janus.openmm_wrapper.keep_atoms(mod_int, qm_atm)
    janus.openmm_wrapper.keep_atoms(mod_str, ['N'])
    atom_int = mod_int.topology.getNumAtoms()
    atom_str = mod_str.topology.getNumAtoms()
    assert atom_int == len(qm_atm)
    assert atom_str == 3

@pytest.mark.datafiles('tests/examples/test_openmm/ala_water.pdb')
def test_delete_residues(datafiles):
    """
    Function to test delete_residues function.
    """
    sys, pdb = create_system(datafiles, 'ala_water.pdb')
    mod_str = janus.openmm_wrapper.create_openmm_modeller(pdb)
    mod_int = janus.openmm_wrapper.create_openmm_modeller(pdb)
    janus.openmm_wrapper.delete_residues(mod_str, ['ALA'])
    janus.openmm_wrapper.delete_residues(mod_int, [0,1,2])
    res_str = mod_str.topology.getNumResidues()
    res_int = mod_int.topology.getNumResidues()
    assert res_str == 28 
    assert res_int == 28 


@pytest.mark.datafiles('tests/examples/test_openmm/ala_water.pdb')
def test_delete_atoms(datafiles):
    """
    Function to test delete_atoms function.
    """
    qm_atm = [0, 1, 2, 3, 4, 5, 6]
    sys, pdb = create_system(datafiles, 'ala_water.pdb')
    mod_str = janus.openmm_wrapper.create_openmm_modeller(pdb)
    mod_int = janus.openmm_wrapper.create_openmm_modeller(pdb)
    janus.openmm_wrapper.delete_atoms(mod_int, qm_atm)
    janus.openmm_wrapper.delete_atoms(mod_str, ['N'])
    atom_int = mod_int.topology.getNumAtoms()
    atom_str = mod_str.topology.getNumAtoms()
    assert atom_int == 117 - len(qm_atm)
    assert atom_str == 114


@pytest.mark.datafiles('tests/examples/test_openmm/ala_water.pdb')
def test_delete_water(datafiles):
    """
    Function to test delete_water function.
    """
    sys, pdb = create_system(datafiles, 'ala_water.pdb')
    mod = janus.openmm_wrapper.create_openmm_modeller(pdb)
    janus.openmm_wrapper.delete_water(mod)
    res = mod.topology.getNumResidues()
    atom = mod.topology.getNumAtoms()
    assert res == 3
    assert atom == 33 
