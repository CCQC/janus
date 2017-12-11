"""
This is the driver module
"""
from . import parser
from .system import System

def energy():
    
    # Initialize system
    sys = parser.parse_input(input_file)
    sys.get_qmmm_energy()
