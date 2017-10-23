import psi4
"""
This module is a wrapper that calls Psi4 to obtain QM information
"""

def get_psi4_energy(sys):
    """
    Call Psi4 to obtain the energy and wavefunction of the QM region 
    and store it in a system object.

    Parameters
    ----------
    sys : System class object

    Returns
    -------
    none
    
    Examples 
    --------
    get_psi4_energy(system)
    """

    mol = psi4.geometry(sys.qm_molecule)
    psi4.set_options(sys.qm_param)
    energy, wavefunction = psi4.energy(sys.qm_method, return_wfn=True)
    sys.qm_energy = energy
