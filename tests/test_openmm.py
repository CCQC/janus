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


@pytest.mark.datafiles('tests/examples/test_openmm/input_protein.pdb')
def test_get_openmm_energy(datafiles):
    """
    Function to test that get_openmm_energy is getting energy correctly.
    """
    sys, pdb = create_system(datafiles, 'input_protein.pdb')
    sim = janus.openmm_wrapper.create_openmm_simulation(sys, 
                                                        pdb.topology,
                                                        pdb.positions)
    potential, kinetic = janus.openmm_wrapper.get_state_info(sim)
    energy = potential + kinetic
    assert np.allclose(energy._value, -495.6537120883586)


@pytest.mark.datafiles('tests/examples/test_openmm/input.pdb')
def test_keep_residues(datafiles):
    """
    Function to test keep_residue function.
    """
    sys, pdb = create_system(datafiles, 'input.pdb')
    mod = janus.openmm_wrapper.create_openmm_modeller(pdb)
    janus.openmm_wrapper.keep_residues(mod, ['HOH'])
    res = mod.topology.getNumResidues()
    assert res == 2761

@pytest.mark.datafiles('tests/examples/test_openmm/input.pdb')
def test_keep_atoms(datafiles):
    """
    Function to test keep_residue function.
    """
    qm_atm = [0, 1, 2, 3, 4, 5, 6, 7]
    sys, pdb = create_system(datafiles, 'input.pdb')
    mod = janus.openmm_wrapper.create_openmm_modeller(pdb)
    janus.openmm_wrapper.keep_atoms(mod, qm_atm)
    atom = mod.topology.getNumAtoms()
    assert atom == len(qm_atm)

@pytest.mark.datafiles('tests/examples/test_openmm/input.pdb')
def test_delete_residues(datafiles):
    """
    Function to test delete_qm_residues function.
    """
    qm_res = [0, 1, 2, 3, 4, 5, 6, 7]
    sys, pdb = create_system(datafiles, 'input.pdb')
    mod = janus.openmm_wrapper.create_openmm_modeller(pdb)
    janus.openmm_wrapper.delete_residues(mod, qm_res)
    res = mod.topology.getNumResidues()
    assert res == 2798 - len(qm_res)


@pytest.mark.datafiles('tests/examples/test_openmm/input.pdb')
def test_delete_atoms(datafiles):
    """
    Function to test delete_qm_atoms function.
    """
    qm_atm = [0, 1, 2, 3, 4, 5, 6, 7]
    sys, pdb = create_system(datafiles, 'input.pdb')
    mod = janus.openmm_wrapper.create_openmm_modeller(pdb)
    janus.openmm_wrapper.delete_atoms(mod, qm_atm)
    atom = mod.topology.getNumAtoms()
    assert atom == 8867 - len(qm_atm)


@pytest.mark.datafiles('tests/examples/test_openmm/input.pdb')
def test_delete_water(datafiles):
    """
    Function to test delete_water function.
    """
    sys, pdb = create_system(datafiles, 'input.pdb')
    mod = janus.openmm_wrapper.create_openmm_modeller(pdb)
    janus.openmm_wrapper.delete_water(mod)
    res = mod.topology.getNumResidues()
    atom = mod.topology.getNumAtoms()
    assert res == 37
    assert atom == 584
