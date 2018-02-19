"""
This is the driver module
"""
from . import parser
from .system import System
from .qm_wrapper import QM_wrapper 
from .mm_wrapper import MM_wrapper 

with open('input.json') as parameter_file:
    parameters = json.load(parameter_file)

system = System(parameters['qmmm'], parameters['qm'], parameters['mm'])

def additive():
    """
    Gets energies of needed components and computes
    a qm/mm energy with a specified embedding method using
    an additive scheme
    """

    # Get MM energy on MM region
    mm_wrapper = MM_wrapper(system, "OpenMM")
    system._ss = mm_wrapper.get_ss()

    # Get nonbonded MM energy on PS-SS interaction
    system._ps_ss = mm_wrapper.get_ps_ss()

    # Get QM energy
    qm_wrapper = QM_wrapper(system, "Psi4")
    system._qm = qm_wrapper.get_qm()

    # Compute total QM/MM energy based on additive scheme
    system.qmmm_energy = system._ss['energy']\
                        + system._ps_ss['energy']\
                        + system._qm['energy']

def subtractive():
    """
    Gets energies of needed components and computes
    a qm/mm energy with a subtractive mechanical embedding scheme
    """

    # Get MM energy on whole system
    mm_wrapper = MM_wrapper(system, "OpenMM")
    system._es = mm_wrapper.get_es()

    # Get MM energy on QM region
    system._ps = mm_wrapper.get_ps()

    # Get QM energy
    qm_wrapper = QM_wrapper(system, "Psi4")
    system._qm = qm_wrapper.get_qm()

    # Compute the total QM/MM energy based on
    # subtractive Mechanical embedding
    system.qmmm_energy = system._es['energy']\
                        - system._ps['energy']\
                        + system._qm['energy']
