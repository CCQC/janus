"""
This is the qmmm driver module
"""
import json
from . import parser
from .system import System
from .qm_wrapper import QM_wrapper 
from .psi4_wrapper import Psi4_wrapper 
from .mm_wrapper import MM_wrapper 
from .openmm_wrapper import OpenMM_wrapper 
from .qmmm import QMMM


def load_system(filename):

    with open(filename) as parameter_file:
        parameters = json.load(parameter_file)

    system = System(parameters['qmmm'], parameters['qm'], parameters['mm'])
    
    return system

def initialize_wrappers(system):
    """
    Initializes the programs to use for computations
    """

    # create qm_wrapper object
    if system.qm_program == "Psi4":
        qm_wrapper = Psi4_wrapper(system)
    else:
    # add other options for qm program here
        pass

    # create mm_wrapper object
    if system.mm_program == "OpenMM":
        mm_wrapper = OpenMM_wrapper(system)
    else:
    # add other options for mm program here
        pass
    
    return mm_wrapper, qm_wrapper

def MD_time_step():
    
    mm_wrapper, qm_wrapper = initialize_wrappers(system)
    
    mm_wrapper.initialize_system()

    for step in range(system.steps):

        mm_wrapper.integration(get_forces=True)
        forces = QMMM.compute_force(scheme=additive)
        mm_wrapper.integration(take_forces=True, take_step=True, force=forces)
        
        
        
        
        

        








