"""
Testing for the openmm_wrapper module
"""
import pytest
import janus
import numpy as np
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *


def create_system(datafiles, filename):
    """
    Function for creating a Janus system with
    mm options used for testing.
    Tests on input_protein.pdb found in tests/examples/test_openmm
    """
    sys = janus.system.System()
    path = str(datafiles)
    sys.mm_pdb_file = os.path.join(path, filename)
    sys.mm_forcefield = 'amber99sb.xml'
    sys.mm_forcefield_water = 'tip3p.xml'
    sys.mm_nonbond_method = PME
    sys.mm_nonbond_cutoff = 1*nanometer
    sys.mm_constraints = HBonds
    sys.mm_temp = 300*kelvin
    sys.qm_res = [0, 1, 2, 3, 4, 5, 6, 7]
    sys.qm_atm = [0, 1, 2, 3, 4, 5, 6, 7]

    janus.openmm_wrapper.create_openmm_system(sys)
    return sys


@pytest.mark.datafiles('tests/examples/test_openmm/input_protein.pdb')
def test_get_openmm_energy(datafiles):
    """
    Function to test that get_openmm_energy is getting energy correctly.
    """
    sys = create_system(datafiles, 'input_protein.pdb')
    janus.openmm_wrapper.create_openmm_simulation(sys)
    janus.openmm_wrapper.get_openmm_energy(sys)
    energy = sys.mm_potential_e + sys.mm_kinetic_e
    assert np.allclose(energy._value, -495.6537120883586)


@pytest.mark.datafiles('tests/examples/test_openmm/input.pdb')
def test_keep_residue(datafiles):
    """
    Function to test keep_residue function.
    """
    sys = create_system(datafiles, 'input.pdb')
    mod = janus.openmm_wrapper.create_openmm_modeller(sys)
    janus.openmm_wrapper.keep_residue(mod, 'HOH')
    res = mod.topology.getNumResidues()
    assert res == 2761


@pytest.mark.datafiles('tests/examples/test_openmm/input.pdb')
def test_delete_qm_residues(datafiles):
    """
    Function to test delete_qm_residues function.
    """
    sys = create_system(datafiles, 'input.pdb')
    mod = janus.openmm_wrapper.create_openmm_modeller(sys)
    janus.openmm_wrapper.delete_qm_residues(mod, sys.qm_res)
    res = mod.topology.getNumResidues()
    assert res == 2798 - len(sys.qm_res)


@pytest.mark.datafiles('tests/examples/test_openmm/input.pdb')
def test_delete_qm_atoms(datafiles):
    """
    Function to test delete_qm_atoms function.
    """
    sys = create_system(datafiles, 'input.pdb')
    mod = janus.openmm_wrapper.create_openmm_modeller(sys)
    janus.openmm_wrapper.delete_qm_atoms(mod, sys.qm_atm)
    atom = mod.topology.getNumAtoms()
    assert atom == 8867 - len(sys.qm_atm)


@pytest.mark.datafiles('tests/examples/test_openmm/input.pdb')
def test_delete_water(datafiles):
    """
    Function to test delete_water function.
    """
    sys = create_system(datafiles, 'input.pdb')
    mod = janus.openmm_wrapper.create_openmm_modeller(sys)
    janus.openmm_wrapper.delete_water(mod)
    res = mod.topology.getNumResidues()
    atom = mod.topology.getNumAtoms()
    assert res == 37
    assert atom == 584
