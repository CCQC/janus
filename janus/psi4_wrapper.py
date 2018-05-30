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

    def qm_info(self):
        if self._energy is None:
            self.get_energy() 
        if self._gradient is None:
            self.get_gradient()
        self._qm['energy'] = self._energy
        self._qm['gradient'] = self._gradient

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
        psi4.core.clean_options()
        self.set_up_psi4()
        self._energy, self._wavefunction = psi4.energy(self._system.qm_method,
                                                        return_wfn=True)

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
        psi4.core.clean_options()
        self.set_up_psi4()
        G = psi4.gradient(self._system.qm_method)
        self._gradient = np.asarray(G)

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
        sys = self._system
        # psi4.core.set_output_file('output.dat', True)
        psi4.core.be_quiet()

        if 'no_reorient' not in sys.qm_positions:
            sys.qm_positions += 'no_reorient \n '
        if 'no_com' not in sys.qm_positions:
            sys.qm_positions += 'no_com \n '

        # make sure this is in nm
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
            # check to make sure positions are in angstroms
        charges = []

        if link is True:
            ss = self._system.second_subsys
            for i, chrg in enumerate(ss['charges']): 
                charges.append([chrg, ss['positions'][i][0], ss['positions'][i][1], ss['positions'][i][2]])
        
        else:

            es = self._system.entire_sys
            mm_index = self._system.boundary_info['mm_index']
            bonds = self._system.boundary_info['bonds_to_mm']
            
            # get q0
            q0 = es['charges'][mm_index] / len(bonds)
            print(q0)
            print(es['charges'][mm_index])
            print(len(bonds))
            print(bonds)
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
        
        pos = []
    
        for bond in bonds:
            new_pos = (positions[bond] + positions[mm]) / 2
            pos.append(new_pos)
        
        return pos
            
                
    def get_scf_charges(self):
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
            self._charges = np.asarray(self._wavefunction.atomic_point_charges())
            self._qm['charges'] = self._charges 


    def get_energy_and_charges(self):
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
        psi4.core.clean_options()
        self.set_up_psi4()
        self._energy, self._wavefunction = psi4.prop(self._system.qm_method,
                                properties=[self._system.qm_charge_method],
                                return_wfn=True)
        self._charges = np.asarray(self._wavefunction.atomic_point_charges())
        self._qm['charges'] = self._charges 


