import psi4
import numpy as np
from . import qm_wrapper
"""
This module is a wrapper that calls Psi4 to obtain QM information
"""

class Psi4_wrapper(QM_wrapper):

    def __init__():
        Super().__init__(system, "Psi4")
        self._qm['energy'] = None
        self._qm['gradient'] = None

def qm_info(self):
    if self._qm['energy'] is None:
        self.get_energy() 
    if self._qm['gradient'] is None:
        self.get_gradient()

def get_energy(self):
    """
    Calls Psi4 to obtain the energy  of the QM region

    Parameters
    ----------
    molecule : a string of molecule parameters in xyz
    param : a dictionary of psi4 parameters
    method : a string of the desired QM method

    Returns
    -------
    An energy

    Examples
    --------
    E = get_psi4_energy(mol, qm_param, 'scf')
    """
    psi4.core.clean()
    set_up_psi4()
    self._energy, self._wavefunction = psi4.energy(self._system.qm_method,
                                       return_wfn=True)
    self._qm['energy'] = self._energy


def get_gradient(self):
    """
    Calls Psi4 to obtain the energy  of the QM region
    and saves it as a numpy array to the passed
    system object as system.qm_gradient

    Parameters
    ----------
    system : a system object containing molecule,
             method, and parameter information

    Returns
    -------
    None

    Examples
    --------
    get_psi4_gradient(system)
    """
    psi4.core.clean()
    set_up_psi4()
    # G, wfn = psi4.gradient(method, return_wfn=True)
    G = psi4.gradient(self._system.qm_method)
    self._qm_gradient = np.asarray(G)
    self._qm['gradient'] = self._qm_gradient

def set_up_psi4(self):
    """
    Sets up a psi4 computation

    Parameters
    ----------
    molecule : a str of molecule parameters
    parameters : A dictionary of psi4 parameters

    Returns
    -------
    None

    Examples
    --------
    set_up_psi4(sys.molecule, sys.parameters)
    """
    ss = self._system._ss
    # psi4.core.set_output_file('output.dat', True)
    psi4.core.be_quiet()

    self._system.qm_positions += 'no_reorient \n '
    self._system.qm_positions += 'no_com \n '
    mol = psi4.geometry(self._system.qm_positions)

    psi4.set_options(self._system.qm_parameters)

    if self._system.embedding_method=='Electrostatic':
        Chrgfield = psi4.QMMM()
        for i in range(len(ss['charge'])):
            Chrgfield.extern.addCharge(ss['charge'][i], ss['positions'][i][0], ss['positions'][i][1], ss['positions'][i][2])
        psi4.core.set_global_option_python('EXTERN', Chrgfield.extern)
            
def get_scf_charge(self):
    """
    Calls Psi4 to obtain the charges on each atom given
    and saves it as a numpy array to the passed system
    object as system.qm_charges.
    This method works well for SCF wavefunctions. For
    correlated levels of theory, it is advised that
    get_psi4_properties() be used instead.

    Parameters
    ----------
    system : a system object containing molecule,
             method, and parameter information

    Returns
    -------
    None

    Examples
    --------
    get_psi4_charge(system)
    """
    if self._wavefunction is not None:
        psi4.oeprop(self._wavefunction, self._system.qm_charge_method)
    self._qm_charges = np.asarray(self._wavefunction.atomic_point_charges())
    self._qm['charges'] = self._qm_charges 


def get_energy_with_charge(self):
    """
    Calls Psi4 to obtain the charges on each atom using the
    property() function and saves it as a numpy array to
    the system as system.qm_charges.

    Parameters
    ----------
    system : a system object containing molecule,
             method, and parameter information

    Returns
    -------
    None

    Examples
    --------
    get_psi4_properties(system)
    """
    psi4.core.clean()
    set_up_psi4()
    self._energy, self._wavefunction = psi4.prop(self._system.qm_method,
                            properties=[self._system.qm_charge_method],
                            return_wfn=True)
    self._qm_charges = np.asarray(self._wavefunction.atomic_point_charges())
    self._qm['charges'] = self._qm_charges 


