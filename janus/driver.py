"""
This is the driver module
"""
import json
from . import parser
from .system import System
from .qm_wrapper import QM_wrapper 
from .psi4_wrapper import Psi4_wrapper 
from .mm_wrapper import MM_wrapper 
from .openmm_wrapper import OpenMM_wrapper 


def load_system(filename):

    with open(filename) as parameter_file:
        parameters = json.load(parameter_file)

    system = System(parameters['qmmm'], parameters['qm'], parameters['mm'])

    return system

def additive(system):
    """
    Gets energientire_sys of needed components and computentire_sys
    a qm/mm energy with a specified embedding method using
    an additive scheme
    """

    #need to add if these things are none then do the following

    # Get MM energy on MM region
    if not system.second_subsys:
        mm_wrapper = OpenMM_wrapper(system)
        system.second_subsys = mm_wrapper.get_second_subsys()

    # Get nonbonded MM energy on PS-SS interaction
    if not system.boundary:
        system.boundary = mm_wrapper.get_boundary()

    # Get QM energy
    if not system.qm:
        # get QM positions from pdb
        if system.qm_positions == None:
            system.qm_positions = mm_wrapper.get_qm_positions() 
        qm_wrapper = Psi4_wrapper(system)
        system.qm = qm_wrapper.get_qm()

    # Compute total QM/MM energy based on additive scheme
    system.qmmm_energy = system.second_subsys['energy']\
                        + system.boundary['energy']\
                        + system.qm['energy']

def subtractive(system):
    """
    Gets energientire_sys of needed components and computentire_sys
    a qm/mm energy with a subtractive mechanical embedding scheme
    """

    # Get MM energy on whole system
    if not system.entire_sys:
        mm_wrapper = OpenMM_wrapper(system)
        system.entire_sys = mm_wrapper.get_entire_sys()

    # Get MM energy on QM region
    if not system.primary_subsys:
        system.primary_subsys = mm_wrapper.get_primary_subsys()

    # Get QM energy
    if not system.qm:
        # get QM positions from pdb
        if system.qm_positions == None:
            system.qm_positions = mm_wrapper.get_qm_positions() 
        qm_wrapper = Psi4_wrapper(system)
        system.qm = qm_wrapper.get_qm()

    # Compute the total QM/MM energy based on
    # subtractive Mechanical embedding
    system.qmmm_energy = system.entire_sys['energy']\
                        - system.primary_subsys['energy']\
                        + system.qm['energy']
