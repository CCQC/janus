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
from .aqmmm import AQMMM
from .oniom_xs import ONIOM_XS


def load_system(filename):

    with open(filename) as parameter_file:
        config = json.load(parameter_file)
    
    return config


def initialize_wrappers(config):
    """
    Initializes the programs to use for computations
    """

    # create qm_wrapper object
    if config['qm_program'] == "Psi4":
        qm_wrapper = Psi4_wrapper(config)
    else:
    # add other options for qm program here
        print("Only Psi4 currently available")

    # create mm_wrapper object
    if config['mm_program'] == "OpenMM":
        mm_wrapper = OpenMM_wrapper(config)
    else:
    # add other options for mm program here
        print("Only OpenMM currently available")


    if config['aqmmm_scheme'] is None:
        qmmm = QMMM(config, qm_wrapper, mm_wrapper)
    elif config['aqmmm_scheme'] == 'ONIOM-XS'
        qmmm = ONIOM_XS(config, qm_wrapper, mm_wrapper)
    else:
        print("Only ONIOM_XS currently implemented")

    # initialize mm_wrapper with information about initial system
    mm_wrapper.initialize()
    
    return mm_wrapper, qmmm

def run_janus(config):

    # initialize wrappers
    mm_wrapper, qmmm = initialize_wrappers(config)


    for step in range(config['steps']):

        #get MM information for entire system
        main_info = mm_wrapper.get_main_info()

        qmmm.run_qmmm(main_info)
        
        # get aqmmm forces 
        forces = qmmm.get_forces()
    
        # feed forces into md simulation and take a step
        # make sure positions are updated so that when i get information on entire system 
        # getting it on the correct one
        mm_wrapper.take_step(force=forces)

