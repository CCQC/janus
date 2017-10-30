"""
This is the driver module
"""
from . import psi4_wrapper as pw
from . import openmm_wrapper as ow
from . import parser
from .system import System


def initialize_system(input_file):
    """
    This function initializes the system with a given input file

    Parameters
    ----------
    input_file : a text file with input parameters

    Returns
    -------
    A System class

    Examples
    --------
    sys = initialize_system(input_file)

    ***maybe put in systems class?
    """
    sys = parser.parse_input(input_file)
    return sys


def qm_energy(sys):
    """
    This function takes a system and calls Psi4 to get the QM energy.
    ***not sure if this is the best approach. Is system passing okay?
    """
    sys.qm_energy = pw.get_psi4_energy(sys.qm_molecule,
                                       sys.qm_method, sys.qm_param)


def mm_energy(sys):
    """
    This function takes a system and calls OpenMM to get MM energy.
    Need to figure out how to use default when no user defined parameters
    """
    mm_system, pdb = ow.create_openmm_system(sys.pdb_file)
    sim = ow.create_openmm_simulation(mm_system, pdb)
    sys.mm_energy = ow.get_openmm_energy(sim)
