import psi4
import numpy as np
from .qm_wrapper import QM_wrapper
"""
This module is a wrapper that calls Psi4 to obtain QM information
"""

class Psi4_wrapper(QM_wrapper):

    def __init__(self, config):

        super().__init__(system, "Psi4")
        self.energy = None
        self.wavefunction = None
        self.gradient = None

    def compute_qm(self):

        self.compute_energy() 
        self.compute_gradient()

    def compute_energy(self):
        """
        Calls Psi4 to obtain the energy and Psi4 wavefunction object of the QM region

        Parameters
        ----------
        None

        Returns
        -------
        Energy, wavefunction

        Examples
        --------
        E = get_psi4_energy()
        """
        psi4.core.clean()
        psi4.core.clean_options()
        self.set_up_psi4()
        self.energy, self.wavefunction = psi4.energy(self.method,
                                                       return_wfn=True)

    def compute_gradient(self):
        """
        Calls Psi4 to obtain the energy  of the QM region
        and saves it as a numpy array self._gradient

        Parameters
        ----------
        None

        Returns
        -------
        None

        Examples
        --------
        get_gradient()
        """
        psi4.core.clean()
        psi4.core.clean_options()
        self.set_up_psi4()
        G = psi4.gradient(self.method)
        self.gradient = np.asarray(G)

    def set_up_psi4(self):
        """
        Sets up a psi4 computation

        Parameters
        ----------
        None

        Returns
        -------
        None

        Examples
        --------
        set_up_psi4()
        """
        # psi4.core.set_output_file('output.dat', True)
        
        # Supress print out
        psi4.core.be_quiet()

        if 'no_reorient' not in self.qm_geometry:
            self.qm_geometry += 'no_reorient \n '
        if 'no_com' not in self.qm_geometry:
            self.qm_geometry += 'no_com \n '

        # make sure this is in angstroms
        mol = psi4.geometry(self.qm_geometry)

        psi4.set_options(self.qm_param)

        if self.charges:
            Chrgfield = psi4.QMMM()
            for charge in self.charges:
                Chrgfield.extern.addCharge(charge[0], charge[1], charge[2], charge[3])
            psi4.core.set_global_option_python('EXTERN', Chrgfield.extern)

            
    def compute_scf_charges(self):
        """
        Calls Psi4 to obtain the charges on each atom given and saves it as a numpy array.
        This method works well for SCF wavefunctions. For correlated levels of theory (e.g., MP2),
        it is advised that get_psi4_properties() be used instead.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Examples
        --------
        get_scf_charge()
        """
        if self.wavefunction is not None:
            psi4.oeprop(self._wavefunction, self.charge_method)
            self.charges = np.asarray(self.wavefunction.atomic_point_charges())
            self.charges = self.charges 


    def compute_energy_and_charges(self):
        """
        Calls Psi4 to obtain the energy, wavefunction, and charges on each atom.
        This method for correlated methods.
        Note: think about passing in wavefunction instead of calling for energy and wavefunction

        Parameters
        ----------
        None

        Returns
        -------
        None

        Examples
        --------
        get_energy_and_charges(system)
        """
        psi4.core.clean()
        psi4.core.clean_options()
        self.set_up_psi4()
        self.energy, self.wavefunction = psi4.prop(self.qm_method,
                                properties=[self.charge_method],
                                return_wfn=True)
        self.charges = np.asarray(self.wavefunction.atomic_point_charges())


    def build_qm_param(self):
        '''
        Builds a dictionary of QM parmeters from input options
        '''
        qm_param = {}
        qm_param['basis'] = self.basis_set
        qm_param['scf_type'] = self.scf_type
        qm_param['guess'] = self.guess_orbitals
        qm_param['reference'] = self.reference
        qm_param['e_convergence'] = self.e_convergence
        qm_param['d_convergence'] = self.d_convergence
        
        self.qm_param = qm_param
