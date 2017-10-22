from . import psi4_wrapper
from . import parser
from .system import System

def qm_energy():
    system = System()
    parser.parse_input('input.dat',system)
    psi4_wrapper.get_psi4_energy(system)
