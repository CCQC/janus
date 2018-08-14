import psi4
import numpy as np
from .qm_wrapper import QM_wrapper
"""
This module is a wrapper that calls Psi4 to obtain QM information
"""

class Psi4_wrapper(QM_wrapper):

    def __init__(self, system):

        super().__init__(system, "Psi4")
        self._energy = None
        self._wavefunction = None
        self._gradient = None

    def qm_info(self, qm_positions):
        if qm_positions:
            self._system.qm_positions = qm_positions

        self.get_energy() 
        self.get_gradient()
        self._qm['energy'] = self._energy
        self._qm['gradients'] = self._gradient

    def get_energy(self):
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
        self._energy, self._wavefunction = psi4.energy(self._system.qm_method,
                                                       return_wfn=True)

    def get_gradient(self):
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
        G = psi4.gradient(self._system.qm_method)
        self._gradient = np.asarray(G)

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
        sys = self._system
        # psi4.core.set_output_file('output.dat', True)
        
        # Supress print out
        psi4.core.be_quiet()

        if 'no_reorient' not in sys.qm_positions:
            sys.qm_positions += 'no_reorient \n '
        if 'no_com' not in sys.qm_positions:
            sys.qm_positions += 'no_com \n '

        # make sure this is in angstroms
        mol = psi4.geometry(sys.qm_positions)

        psi4.set_options(sys.qm_param)

        if sys.embedding_method=='Electrostatic':
            if sys.boundary_treatment == 'link_atom' or None: 
                charges = self.get_external_charges(link=True)
            if sys.boundary_treatment == 'RC': 
                charges = self.get_external_charges(RC=True)
            if sys.boundary_treatment == 'RCD': 
                charges = self.get_external_charges(RCD=True)
            Chrgfield = psi4.QMMM()
            for charge in charges:
                Chrgfield.extern.addCharge(charge[0], charge[1], charge[2], charge[3])
            psi4.core.set_global_option_python('EXTERN', Chrgfield.extern)


    def get_external_charges(self, link=False, RC=False, RCD=False):
        """
        Gets the point charges of atoms from secondary subsystem for electrostatic embedding 

        Note: check to make sure positions are in angstroms

        Parameters
        ----------
        link: a bool specifying whether the link atom scheme is used for determining point charges. 
              default is False.
        RC: a bool specifying whether the RC scheme is used for determining point charges. 
              default is False.
        RCD: a bool specifying whether the RCD scheme is used for determining point charges. 
              default is False.

        Returns
        -------
        A list of charges and cooresponding positions as  xyz coordinates

        Examples
        --------
        get_external_charge(link=True)
        get_external_charge(RC=True)
        """
        charges = []

        if link is True:
            ss = self._system.second_subsys
            for i, chrg in enumerate(ss['charges']): 
                charges.append([chrg, ss['positions'][i][0], ss['positions'][i][1], ss['positions'][i][2]])
        
        # This is for the RC and RCD schemes
        else:

            es = self._system.entire_sys
            mm_index = self._system.boundary_info['mm_index']
            bonds = self._system.boundary_info['bonds_to_mm']
            
            # get q0
            q0 = es['charges'][mm_index] / len(bonds)

            # get positions
            positions = self.get_redistributed_positions(es['positions'], bonds, mm_index)

            if RC is True:
                for i, chrg in enumerate(es['charges']):
                    # add every atom not in qm system or the M1 atom 
                    if i not in self._system.qm_atoms and i != mm_index:
                            charges.append([chrg, es['positions'][i][0], es['positions'][i][1], es['positions'][i][2]])
                for pos in positions:
                    charges.append([q0, pos[0], pos[1], pos[2]])

            elif RCD is True:
                q0_RCD = q0 * 2

                for i, chrg in enumerate(es['charges']):
                    # add every atom not in qm system or the M1 atom 
                    if i not in self._system.qm_atoms and i != mm_index:
                        if i in bonds:
                        # modified M2 charges in RCD scheme
                            charges.append([chrg - q0, es['positions'][i][0], es['positions'][i][1], es['positions'][i][2]])
                        else:
                            charges.append([chrg, es['positions'][i][0], es['positions'][i][1], es['positions'][i][2]])

                for pos in positions:
                    charges.append([q0_RCD, pos[0], pos[1], pos[2]])
                
        return charges

    def get_redistributed_positions(self, positions, bonds, mm):
        """
        Gets the positions for the redistributed point charges in the RC and RCD schemes

        Parameters
        ----------
        positions: a list of the positions
        bonds: a list of indices of all atoms (in secondary subsystem) bonded to M1  
        mm: the index of M1

        Returns
        -------
        List of positions for the redistributed charges

        Examples
        --------
        get_redistributed_positions(positions=pos, bonds=bond, mm=mm_index)
        """
        
        pos = []
    
        for bond in bonds:
            new_pos = (positions[bond] + positions[mm]) / 2
            pos.append(new_pos)
        
        return pos
            
                
    def get_scf_charges(self):
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
        if self._wavefunction is not None:
            psi4.oeprop(self._wavefunction, self._system.qm_charge_method)
            self._charges = np.asarray(self._wavefunction.atomic_point_charges())
            self._qm['charges'] = self._charges 


    def get_energy_and_charges(self):
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
        self._energy, self._wavefunction = psi4.prop(self._system.qm_method,
                                properties=[self._system.qm_charge_method],
                                return_wfn=True)
        self._charges = np.asarray(self._wavefunction.atomic_point_charges())
        self._qm['charges'] = self._charges 


