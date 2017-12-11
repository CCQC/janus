import psi4
import numpy as np
"""
This module is a wrapper that calls Psi4 to obtain QM information
"""


def get_psi4_energy(molecule, param, method):
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
    set_up_psi4(molecule, param)
    energy = psi4.energy(method)
   # energy, wavefunction = psi4.energy(method,
   #                                    return_wfn=True)
    return energy


def get_psi4_charge(system):
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
    get_psi4_wavefunction(system)
    psi4.oeprop(system.qm_wfn, system.qm_charge_method)
    system.qm_charges = np.asarray(system.qm_wfn.atomic_point_charges())


def get_psi4_properties(system):
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
    set_up_psi4(system.qm_molecule, system.qm_param)
    energy, wfn = psi4.prop(system.qm_method,
                            properties=[system.qm_charge_method],
                            return_wfn=True)
    system.qm_charges = np.asarray(wfn.atomic_point_charges())


def get_psi4_gradient(system):
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
    set_up_psi4(system.qm_molecule, system.qm_param)
    # G, wfn = psi4.gradient(method, return_wfn=True)
    G = psi4.gradient(system.qm_method)
    system.qm_gradient = np.asarray(G)


def get_psi4_qmmm():
    pass


def get_psi4_wavefunction(system):
    """
    Calls Psi4 to obtain the wavefunction of the QM region

    Parameters
    ----------
    system : a system object containing molecule,
             method, and parameter information

    Returns
    -------
    A Psi4 wavefunction object

    Examples
    --------
    get_psi4_wavefunction(system)
    """
    set_up_psi4(system.qm_molecule, system.qm_param)
    energy, wavefunction = psi4.energy(system.qm_method,
                                       return_wfn=True)
    system.qm_wfn = wavefunction


def set_up_psi4(molecule, parameters):
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
    # psi4.core.set_output_file('output.dat', True)
    psi4.core.be_quiet()
    mol = psi4.geometry(molecule)
    psi4.set_options(parameters)
