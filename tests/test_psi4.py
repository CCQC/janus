"""
Testing for psi4_wrapper.py module
"""
import pytest
import janus
import psi4
import numpy as np

def test_get_psi4_energy():
    """
    Function to test get_psi4_energy is getting energy correctly
    """
    
    qm_mol = """
0 1
O
H 1 R
H 1 R 2 A
R = 1.0
A = 104.5
symmetry c1"""

    parameters = \
    {
        'basis' : 'STO-3G',
        'scf_type' : 'df',
        'guess' : 'sad',
        'reference' : 'rhf',
        'e_convergence' : 1e-8,
        'd_convergence' : 1e-8
    }

    method = 'scf'
    
    sys = janus.system.System()
    sys.qm_molecule, sys.qm_param = qm_mol, parameters
    sys.qm_method = method

    janus.psi4_wrapper.get_psi4_energy(sys)
    
    mol = psi4.geometry(qm_mol)
    psi4.set_options(parameters)
    energy = psi4.energy(method, molecule = mol)
    
    assert np.allclose(sys.qm_energy, energy)
