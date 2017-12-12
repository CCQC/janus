"""
Testing for the System class
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
    pdb_file = os.path.join(path, filename)

    qm_mol = """O     0.123   3.593   5.841 
 H    -0.022   2.679   5.599 
 H     0.059   3.601   6.796 
 O     0.017   6.369   7.293 
 H    -0.561   5.928   6.669 
 H     0.695   6.771   6.749 
 no_reorient 
 no_com 
 """

    parameters = \
        {
            'basis': 'STO-3G',
            'scf_type': 'df',
            'guess': 'sad',
            'reference': 'rhf',
            'e_convergence': 1e-8,
            'd_convergence': 1e-8
        }

    sys = janus.system.System(qm_param=parameters, qm_method='scf',
                              qm_molecule=qm_mol, qm_atoms=[0, 1, 2, 3, 4, 5],
                              mm_pdb_file=pdb_file)
    return sys


@pytest.mark.datafiles('tests/examples/test_openmm/water.pdb')
def test_get_openmm_state_info(datafiles):
    """
    Function to test get_openmm_energy function
    of systems class
    """
    sys = create_system(datafiles, 'water.pdb')
    sys.get_openmm_state_info()
    assert np.allclose(sys.mm_Te, -0.010571307078971566)
    assert np.allclose(sys.mm_Ke, 8.414710565572852e-06)
    assert np.allclose(sys.mm_tot_energy, -0.010562892368405992)


@pytest.mark.datafiles('tests/examples/test_openmm/water.pdb')
def test_get_mm_qm_energy(datafiles):
    """
    Function to test get_mm_qm_energy function
    of systems class
    """
    sys = create_system(datafiles, 'water.pdb')
    sys.get_mm_qm_energy()
    assert np.allclose(sys.mod_Te, -0.005239479792864975)
    assert np.allclose(sys.mod_Ke, 3.0368070644980085e-06)
    assert np.allclose(sys.mm_qm_energy, -0.005236442985800477)


@pytest.mark.datafiles('tests/examples/test_openmm/water.pdb')
def test_get_qmmm_energy(datafiles):
    """
    Function to test get_qmmm_energy function
    of systems class given the mm energy and the mm energy
    of the qm region
    """
    sys = create_system(datafiles, 'water.pdb')
    sys.mm_tot_energy = -0.010562892368405992
    sys.mm_qm_energy = -0.005236442985800477
    sys.get_qmmm_energy()
    assert np.allclose(sys.qm_energy, -149.92882700821423)
    assert np.allclose(sys.qmmm_energy, -149.93415345759684)

@pytest.mark.datafiles('tests/examples/test_openmm/water.pdb')
def test_make_qm_molecule(datafiles):
    """
    Function to test get_qmmm_energy function
    of systems class given the mm energy and the mm energy
    of the qm region
    """
    sys = create_system(datafiles, 'water.pdb')
    mol = sys.make_qm_molecule()
    assert mol == sys.qm_molecule
