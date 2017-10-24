"""
This is the driver module 
"""
from . import psi4_wrapper
from . import parser
from .system import System


def qm_energy():
    """
    This function instantiates a system object,
    stores input information with parse_input, and
    calls Psi4 to get the QM energy.
    """
    sys = System()
    parser.parse_input('input.dat', system)
    sys.qm_energy = psi4_wrapper.get_psi4_energy(sys.molecule, sys.qm_method, sys.qm_param)

def mm_energy():
    pass
