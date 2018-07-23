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


def run_adaptive(filename):
   
    # initialize system 
    system = load_system(filename)

    # initialize wrappers
    mm_wrapper, qm_wrapper = initialize_wrappers(system)

    # with openmm wrapper,
    # this creates 2 openmm objects containing entire system
    # one for computing forces and one for time step integration
    mm_wrapper.initialize_system()

    aqmmm = AQMMM()

    for step in range(system.steps):

        # get MM information for entire system
        entire_sys = mm_wrapper.get_entire_sys()

        # get the partitions for each qmmm computation
        paritions = aqmmm.partition(entire_sys, system.aqmmm_partition_scheme)
        for partition in partitions:
            qmmm = qmmm.get_info(scheme=system.qmmm_scheme, entire_sys, partition)
            aqmmm.save(partition.ID, qmmm.forces, qmmm.energy)
            
        # get aqmmm forces 
        forces = adaptive.get_info(scheme=system.aqmmm_scheme)
    
        # feed forces into md simulation and take a step
        # make sure positions are updated so that when i get information on entire system 
        # getting it on the correct one
        mm.wrapper.take_step(force=forces)

def run_qmmm():
# have this as part of run_adaptive?

    # initialize system 
    system = load_system(filename)

    # initialize wrappers
    mm_wrapper, qm_wrapper = initialize_wrappers(system)

    # with openmm wrapper,
    # this creates 2 openmm objects containing entire system
    # one for computing forces and one for time step integration
    mm_wrapper.initialize_system()

    for step in range(system.steps):

        # get MM information for entire system
        entire_sys = mm_wrapper.get_entire_sys()

        qmmm = qmmm.get_info(scheme=system.qmmm_scheme, entire_sys)
            
        # get aqmmm forces 
        forces = adaptive.get_info(scheme=system.aqmmm_scheme)
    
        # feed forces into md simulation and take a step
        # make sure positions are updated so that when i get information on entire system 
        # getting it on the correct one
        mm.wrapper.take_step(force=forces)

