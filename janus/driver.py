"""
This is the driver module
"""
from . import parser
from .system import System

def energy(input_file):
    
    # Initialize system
    sys = parser.parse_input(input_file)
    sys.make_qm_molecule()
    sys.get_qmmm_energy()
