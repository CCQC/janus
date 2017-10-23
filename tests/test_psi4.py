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
        'scf_type' : 'pk',
        'reference' : 'rhf',
        'e_convergence' : 1e-8,
        'd_convergence' : 1e-8
    }

    system = janus.system.System(parameters, 'scf', qm_mol)
    janus.psi4_wrapper.get_psi4_energy(system)
    
    mol = psi4.geometry("""
                        0 1
                        O
                        H 1 R
                        H 1 R 2 A
                        R = 1.0
                        A = 104.5
                        symmetry c1
                        """)

    psi4.set_options({'basis' : 'STO-3G',
                      'scf_type' : 'pk',
                      'reference' : 'rhf',
                      'e_convergence' : 1e-8,
                      'd_convergence' : 1e-8})

    energy = psi4.energy('scf', molecule = mol)
    
    assert np.allclose(system.qm_energy, energy)
    
    
