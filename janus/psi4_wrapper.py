import psi4
"""
This module is a wrapper that calls Psi4 to obtain QM information
"""


def get_psi4_energy(molecule, method, parameters):
    """
    Calls Psi4 to obtain the energy and wavefunction of the QM region
    and store it in a system object.

    Parameters
    ----------
    molecule : a str of molecule parameters
    method : a str of desired method
    parameters : A dictionary of psi4 parameters
    sys : System class object

    Returns
    -------
    A psi4 energy

    Examples
    --------
    get_psi4_energy(sys.molecule, sys.qm_method, sys.qm_param)
    """

    mol = psi4.geometry(molecule)
    psi4.set_options(parameters)
    #energy, wavefunction = psi4.energy(sys.qm_method, return_wfn=True)
    energy = psi4.energy(method)
    return energy
