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
    Gets energies of needed components and computes
    a qm/mm energy with a specified embedding method using
    an additive scheme
    """

    # Get MM energy on MM region
    mm_wrapper = OpenMM_wrapper(system)
    system.ss = mm_wrapper.get_ss()

    print(system.ss['energy'])
    # Get nonbonded MM energy on PS-SS interaction
    system.ps_ss = mm_wrapper.get_ps_ss()
    print(system.ps_ss['energy'])

    # get QM positions from pdb
    system.qm_positions = mm_wrapper.get_qm_positions() 
    # Get QM energy
    qm_wrapper = Psi4_wrapper(system)
    system.qm = qm_wrapper.get_qm()

    print(system.qm['energy'])
    # Compute total QM/MM energy based on additive scheme
    system.qmmm_energy = system.ss['energy']\
                        + system.ps_ss['energy']\
                        + system.qm['energy']

def subtractive(system):
    """
    Gets energies of needed components and computes
    a qm/mm energy with a subtractive mechanical embedding scheme
    """

    # Get MM energy on whole system
    mm_wrapper = OpenMM_wrapper(system)
    system.es = mm_wrapper.get_es()

    # Get MM energy on QM region
    system.ps = mm_wrapper.get_ps()

    # get QM positions from pdb
    system.qm_positions = mm_wrapper.get_qm_positions() 

    # Get QM energy
    qm_wrapper = Psi4_wrapper(system)
    system.qm = qm_wrapper.get_qm()

    # Compute the total QM/MM energy based on
    # subtractive Mechanical embedding
    system.qmmm_energy = system.es['energy']\
                        - system.ps['energy']\
                        + system.qm['energy']
