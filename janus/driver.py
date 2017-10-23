from . import psi4_wrapper
from . import parser
from .system import System

"""
This module is the driver 
"""

def qm_energy():
    """
    This function instantiates a system object, 
    stores input information with parse_input, and 
    calls Psi4 to get the QM energy. 
    """
    system = System()
    parser.parse_input('input.dat',system)
    psi4_wrapper.get_psi4_energy(system)
