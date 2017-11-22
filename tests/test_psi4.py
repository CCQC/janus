"""
Testing for psi4_wrapper.py module
"""
import pytest
import janus
import psi4
import numpy as np


def create_system():
    """
    Function for creating a Janus system with
    qm options used for testing
    """
    qm_mol = """
0 1
O
H 1 R
H 1 R 2 A
R = 1.0
A = 104.5
symmetry c1
no_reorient
no_com"""

    parameters = \
        {
            'basis': 'STO-3G',
            'scf_type': 'df',
            'guess': 'sad',
            'reference': 'rhf',
            'e_convergence': 1e-8,
            'd_convergence': 1e-8
        }

    sys = janus.system.System()
    sys.qm_molecule, sys.qm_param = qm_mol, parameters

    return sys


def test_get_psi4_energy():
    """
    Function to test get_psi4_energy is getting energy correctly
    """
    sys = create_system()

    energy = janus.psi4_wrapper.get_psi4_energy(sys.qm_molecule, sys.qm_param, sys.qm_method)

    assert np.allclose(energy, -74.96475172187111)


def test_get_psi4_charge():
    """
    Function to test get_psi4_charge is getting the charges correctly
    """
    sys = create_system()

    janus.psi4_wrapper.get_psi4_charge(sys)

    charges = np.asarray([-0.32804757, 0.16402379, 0.16402379])

    assert np.allclose(sys.qm_charges, charges)


def test_get_psi4_properties():
    """
    Function to test get_psi4_properties is getting the charges correctly
    """
    sys = create_system()

    janus.psi4_wrapper.get_psi4_properties(sys)

    charges = np.asarray([-0.32804757, 0.16402379, 0.16402379])

    assert np.allclose(sys.qm_charges, charges)


def test_get_psi4_gradient():
    """
    Function to test get_psi4_gradient is getting the gradient correctly
    """
    sys = create_system()

    janus.psi4_wrapper.get_psi4_gradient(sys)

    gradient = np.asarray([[-7.78353232e-16, 5.14539763e-12, 1.77427098e-03],
                           [1.03893379e-15, -1.93698815e-02, -8.87135490e-04],
                           [-2.60580557e-16, 1.93698815e-02, -8.87135492e-04]])

    assert np.allclose(sys.qm_gradient, gradient)
