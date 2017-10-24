"""
Testing for the openmm_wrapper module
"""
import pytest
import janus
import numpy as np
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *


@pytest.mark.datafiles('tests/examples/test_openmm/input_protein.pdb')
def test_get_openmm_energy(datafiles):
    """
    Function to test that get_openmm_energy is getting energy correctly.
    Tests on input_protein.pdb found in tests/examples/test_openmm
    """

    path = str(datafiles)
    sys, pdb = janus.openmm_wrapper.create_openmm_system(os.path.join(path, 'input_protein.pdb'))
    sim = janus.openmm_wrapper.create_openmm_simulation(sys, pdb)
    energy = janus.openmm_wrapper.get_openmm_energy(sim)
    assert np.allclose(energy._value, -495.6537120883586)


@pytest.mark.datafiles('tests/examples/test_openmm/input.pdb')
def test_keep_residue(datafiles):
    """
    Function to test keep_residue function.
    Tests on input.pdb found in tests/examples/test_openmm
    """

    path = str(datafiles)
    pdb = PDBFile(os.path.join(path, 'input.pdb'))
    mod = janus.openmm_wrapper.create_openmm_modeller(pdb)
    janus.openmm_wrapper.keep_residue(mod, 'HOH')
    res = mod.topology.getNumResidues()
    assert res == 2761


@pytest.mark.datafiles('tests/examples/test_openmm/input.pdb')
def test_delete_qm_residues(datafiles):
    """
    Function to test delete_qm_residues function.
    Tests on input.pdb found in tests/examples/test_openmm
    """

    path = str(datafiles)
    pdb = PDBFile(os.path.join(path, 'input.pdb'))
    mod = janus.openmm_wrapper.create_openmm_modeller(pdb)
    qm_res = [0, 1, 2, 3, 4, 5, 6, 7]
    janus.openmm_wrapper.delete_qm_residues(mod, qm_res)
    res = mod.topology.getNumResidues()
    assert res == 2798 - len(qm_res)


@pytest.mark.datafiles('tests/examples/test_openmm/input.pdb')
def test_delete_qm_atoms(datafiles):
    """
    Function to test delete_qm_atoms function.
    Tests on input.pdb found in tests/examples/test_openmm
    """

    path = str(datafiles)
    pdb = PDBFile(os.path.join(path, 'input.pdb'))
    mod = janus.openmm_wrapper.create_openmm_modeller(pdb)
    qm_atm = [0, 1, 2, 3, 4, 5, 6, 7]
    janus.openmm_wrapper.delete_qm_atoms(mod, qm_atm)
    atom = mod.topology.getNumAtoms()
    assert atom == 8867 - len(qm_atm)


@pytest.mark.datafiles('tests/examples/test_openmm/input.pdb')
def test_delete_water(datafiles):
    """
    Function to test delete_water function.
    Tests on input.pdb found in tests/examples/test_openmm
    """

    path = str(datafiles)
    pdb = PDBFile(os.path.join(path, 'input.pdb'))
    mod = janus.openmm_wrapper.create_openmm_modeller(pdb)
    janus.openmm_wrapper.delete_water(mod)
    res = mod.topology.getNumResidues()
    atom = mod.topology.getNumAtoms()
    assert res == 37
    assert atom == 584
