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

        self.get_energy() 
        self.get_gradient()

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

#        if sys.embedding_method=='Electrostatic':
#            if sys.boundary_treatment == 'link_atom' or None: 
#                charges = self.get_external_charges(link=True)
#            if sys.boundary_treatment == 'RC': 
#                charges = self.get_external_charges(RC=True)
#            if sys.boundary_treatment == 'RCD': 
#                charges = self.get_external_charges(RCD=True)
#            Chrgfield = psi4.QMMM()
#            for charge in charges:
#                Chrgfield.extern.addCharge(charge[0], charge[1], charge[2], charge[3])
#            psi4.core.set_global_option_python('EXTERN', Chrgfield.extern)


'''
PUT THIS IN QMMM
'''
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
        # REWRITE!!!!!
#        charges = []
#
#        if link is True:
#            ss = self._system.second_subsys
#            for i, chrg in enumerate(ss['charges']): 
#                charges.append([chrg, ss['positions'][i][0], ss['positions'][i][1], ss['positions'][i][2]])
#        
#        # This is for the RC and RCD schemes
#        else:
#
#            es = self._system.entire_sys
#            mm_index = self._system.boundary_info['mm_index']
#            bonds = self._system.boundary_info['bonds_to_mm']
#            
#            # get q0
#            q0 = es['charges'][mm_index] / len(bonds)
#
#            # get positions
#            positions = self.get_redistributed_positions(es['positions'], bonds, mm_index)
#
#            if RC is True:
#                for i, chrg in enumerate(es['charges']):
#                    # add every atom not in qm system or the M1 atom 
#                    if i not in self._system.qm_atoms and i != mm_index:
#                            charges.append([chrg, es['positions'][i][0], es['positions'][i][1], es['positions'][i][2]])
#                for pos in positions:
#                    charges.append([q0, pos[0], pos[1], pos[2]])
#
#            elif RCD is True:
#                q0_RCD = q0 * 2
#
#                for i, chrg in enumerate(es['charges']):
#                    # add every atom not in qm system or the M1 atom 
#                    if i not in self._system.qm_atoms and i != mm_index:
#                        if i in bonds:
#                        # modified M2 charges in RCD scheme
#                            charges.append([chrg - q0, es['positions'][i][0], es['positions'][i][1], es['positions'][i][2]])
#                        else:
#                            charges.append([chrg, es['positions'][i][0], es['positions'][i][1], es['positions'][i][2]])
#
#                for pos in positions:
#                    charges.append([q0_RCD, pos[0], pos[1], pos[2]])
#                
#        return charges

'''
PUT THE FOLLOWING IN QMMM!!
'''
    #def get_redistributed_positions(self, positions, bonds, mm):
    #    """
    #    Gets the positions for the redistributed point charges in the RC and RCD schemes

    #    Parameters
    #    ----------
    #    positions: a list of the positions
    #    bonds: a list of indices of all atoms (in secondary subsystem) bonded to M1  
    #    mm: the index of M1

    #    Returns
    #    -------
    #    List of positions for the redistributed charges

    #    Examples
    #    --------
    #    get_redistributed_positions(positions=pos, bonds=bond, mm=mm_index)
    #    """
    #    
    #    pos = []
    #
    #    for bond in bonds:
    #        new_pos = (positions[bond] + positions[mm]) / 2
    #        pos.append(new_pos)
    #    
    #    return pos
            
                
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
