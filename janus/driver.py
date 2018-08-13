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

    system = System(parameters['aqmmm'], parameters['qmmm'], parameters['qm'], parameters['mm'])
    
    return system

def initialize_wrappers(system, aqmmm=False):
    """
    Initializes the programs to use for computations
    """

    if aqmmm is False:
        # create qm_wrapper object
        if system.qm_program == "Psi4":
            qm_wrapper = Psi4_wrapper(system)
        else:
        # add other options for qm program here
            print("Only Psi4 currently available")

        # create mm_wrapper object
        if system.mm_program == "OpenMM":
            mm_wrapper = OpenMM_wrapper(system)
        else:
        # add other options for mm program here
            print("Only OpenMM currently available")

        return mm_wrapper, qm_wrapper

    if aqmmm is True:
        if system.aqmmm_scheme == 'ONIOM_XS':
            aqmmm = ONIOM_XS(system.aqmmm_partition_scheme, system.mm_pdb_file)
        else:
            print("Only ONIOM_XS currently implemented")
        
        return aqmmm


def run_adaptive(system):

    # initialize wrappers
    mm_wrapper, qm_wrapper = initialize_wrappers(system)

    aqmm = initialize_wrappers(system, aqmmm=True)

    # with openmm wrapper,
    # this creates 2 openmm objects containing entire system
    # one for computing forces and one for time step integration
    # each individual mm_wrapper gets initial trajectory


    qmmm = QMMM(qm_wrapper, mm_wrapper)

    for step in range(system.steps):

        # get MM information for entire system
        # main_info = mm_wrapper.get_main_info()

        # get the partitions for each qmmm computation
        # the thing passed in will have positions for the trajectory to be updated
        paritions = aqmmm.partition(entire_sys)
        for partition in partitions:
            qmmm.get_info(system.qmmm_scheme, mm_wrapper, partitition=partition)
            aqmmm.save(partition.ID, qmmm.qmmm_forces, qmmm.qmmm_energy)
            
        # get aqmmm forces 
        forces = aqmmm.get_info()
    
        # feed forces into md simulation and take a step
        # make sure positions are updated so that when i get information on entire system 
        # getting it on the correct one
        mm_wrapper.take_step(force=forces)

def run_qmmm(system):
# have this as part of run_adaptive?

    # initialize wrappers
    mm_wrapper, qm_wrapper = initialize_wrappers(system)


    # with openmm wrapper,
    # this creates 2 openmm objects containing entire system
    # one for computing forces and one for time step integration
    trajectory = mm_wrapper.initialize_system()

    qmmm = QMMM(qm_wrapper)

    for step in range(system.steps):

        # get MM information for entire system
        # main_info = mm_wrapper.get_main_info()

        qmmm.get_info(system.qmmm_scheme, mm_wrapper)
            
        # get qmmm forces 
        forces = qmmm.qmmm_forces 
    
        # feed forces into md simulation and take a step
        # make sure positions are updated so that when i get information on entire system 
        # getting it on the correct one
        mm_wrapper.take_step(force=forces)

