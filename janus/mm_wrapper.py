"""
MM wrapper super clasecond_subsys
"""
from abc import ABC, abstractmethod

class MM_wrapper(ABC):

    kjmol_to_au = 1/2625.5002 
    nm_to_angstrom = 10.0000000
    nm_to_bohr = 0.052917999999 
    kj_mol_nm_to_au_bohr = kjmol_to_au*nm_to_bohr 

    def __init__(self, system, program):
        
        self._program = program
        self._system = system
        self._second_subsys = {}
        self._primary_subsys = {}
        self._entire_sys = {}
        self._boundary = {} 
        self._boundary['energy'] = None
        self._qm_positions = None
        self.link_atoms = {}

        super().__init__()

    @abstractmethod
    def second_subsys_info(self):
        pass

    @abstractmethod
    def primary_subsys_info(self, link=False):
        pass

    @abstractmethod
    def entire_sys_info(self):
        pass

    @abstractmethod
    def boundary_info(self):
        pass

    @abstractmethod
    def qm_positions(self):
        pass

    @abtractmethod
    def integration(self, get_forces, take_forces, take_step, force):
        pass


#    @abstractmethod
#    def make_zero_energy(self):
#        pass

    def get_second_subsys(self, coulomb=True):
        """
        Gets the information (energy, position, forces) for the secondary subsystem only.

        Parameters
        ----------
        coulomb: Whether to include coulombic interactions in MM computation.
                 Default is True.

        Returns
        -------
        Dictionary with relevant information
        """

        if not self._second_subsys:
            self.second_subsys_info(coulomb)
        return self._second_subsys

    def get_primary_subsys(self, link=False, coulomb=True):
        """
        Gets the information (energy, position, forces) for the primary subsystem only.
        Note: these additional options might need to be openmm specific, we'll see

        Parameters
        ----------
        link: Whether to include a link atom in the MM computation. 
        coulomb: Whether to include coulombic interactions in MM computation.
                 Default is True.

        Returns
        -------
        Dictionary with relevant information
        """

        if not self._primary_subsys:
            self.primary_subsys_info(link, coulomb)
        return self._primary_subsys
    
    def get_entire_sys(self, coulomb=True):
        """
        Gets the information (energy, position, forces) for the entire system only.

        Parameters
        ----------
        coulomb: Whether to include coulombic interactions in MM computation.
                 Default is True.

        Returns
        -------
        Dictionary with relevant information
        """

        if not self._entire_sys:
            self.entire_sys_info(coulomb)
        return self._entire_sys

    def get_boundary(self, coulomb=True):
        """
        Gets the information for the interaction energy between the primary and secondary subsystem (QM-MM)

        Parameters
        ----------
        coulomb: Whether to include coulombic interactions in MM computation.
                 Default is True.

        Returns
        -------
        Dictionary with relevant information
        """

        if self._boundary['energy'] is None:
            self.boundary_info(coulomb)
        return self._boundary

    def get_qm_positions(self):
        """
        Grabs the positions of the atoms in the primary subsystem. 
        Note: In a MD time step, does the position update or not? Need to make sure this updates

        Parameters
        ----------
        None

        Returns
        -------
        A string with the element and xyz geometry coordinate of atoms of primary subsystem
        """

        if self._qm_positions is None:
            self.qm_positions()
        return self._qm_positions

    def get_boundary_info(self):
        """
        Gets information for any link atoms 

        Parameters
        ----------
        None

        Returns
        -------
        A dictionary with link atom information
        """

        if self.link_atoms:
            return self.link_atoms
