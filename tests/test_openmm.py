import pytest
import janus
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

#FIXTURE_DIR = os.path.join(
#    os.path.dirname(os.path.realpath(__file__)),
#    'tests/examples/test_openmm/',
#    )
#
#ALL_FILES = pytest.mark.datafiles(
#    os.path.join(FIXTURE_DIR, 'input.pdb'),
#    )

def test_get_openmm_energy():
    pass

@pytest.mark.datafiles('tests/examples/test_openmm/input.pdb')
def test_keep_residue(datafiles):
    path = str(datafiles)
    pdb = PDBFile(os.path.join(path, 'input.pdb'))
    mod = janus.openmm_wrapper.create_openmm_modeller(pdb)
    janus.openmm_wrapper.keep_residue(mod, 'HOH')
    res = mod.topology.getNumResidues()
    assert res == 2761

@pytest.mark.datafiles('tests/examples/test_openmm/input.pdb')
def test_delete_qm_residues(datafiles):
    path = str(datafiles)
    pdb = PDBFile(os.path.join(path, 'input.pdb'))
    mod = janus.openmm_wrapper.create_openmm_modeller(pdb)
    qm_res = [0, 1, 2, 3, 4, 5, 6, 7]
    janus.openmm_wrapper.delete_qm_residues(mod, qm_res)
    res = mod.topology.getNumResidues()
    assert res == 2798 - len(qm_res)

@pytest.mark.datafiles('tests/examples/test_openmm/input.pdb')
def test_delete_qm_atoms(datafiles):
    path = str(datafiles)
    pdb = PDBFile(os.path.join(path, 'input.pdb'))
    mod = janus.openmm_wrapper.create_openmm_modeller(pdb)
    qm_atm = [0, 1, 2, 3, 4, 5, 6, 7]
    janus.openmm_wrapper.delete_qm_atoms(mod, qm_atm)
    atom = mod.topology.getNumAtoms()
    assert atom == 8867 - len(qm_atm)

@pytest.mark.datafiles('tests/examples/test_openmm/input.pdb')
def test_delete_water(datafiles):
    path = str(datafiles)
    pdb = PDBFile(os.path.join(path, 'input.pdb'))
    mod = janus.openmm_wrapper.create_openmm_modeller(pdb)
    janus.openmm_wrapper.delete_water(mod)
    res = mod.topology.getNumResidues()
    atom = mod.topology.getNumAtoms()
    assert res == 37
    assert atom == 584
