from abc import ABC, abstractmethod

class MM_wrapper(ABC):

    kjmol_to_au = 1/2625.5002 
    nm_to_angstrom = 10.0000000
    nm_to_bohr = 0.052917999999 
    kjmol_nm_to_au_bohr = kjmol_to_au*nm_to_bohr 
    au_bohr_to_kjmol_nm = 1/kjmol_nm_to_au_bohr

    def __init__(self, param, class_type):
        """
        A super class for all molecular mechanics wrappers

        Parameters
        ----------
        param : dict 
        program : str
            what program to use for QM computations

        """
        self.system_info = param['system_info']
        self.system_info_format = param['system_info_format']
        self.param = param
        self.class_type = class_type
        self.main_info = None
        self.main_charges = None

        super().__init__()

    def get_energy_and_gradient(self, traj, include_coulomb='all', link_atoms=None, minimize=False, charges=None):
        """
        Gets the energy and gradient from a MM computation

        Parameters
        ----------
        traj : MDtraj trajectory object
        include_coulomb : str

            whether to include coulombic interactions. 
            'all' (default) includes coulombic forces for all particles,
            'no_link' excludes coulombic forces for link atoms,
            'only' excludes all other forces for all atoms,
            'none' excludes coulombic forces for all particles.

        link_atoms : list
            indices of link_atoms
        minimize : bool
            whether to return the geometry optimized energy 
        charges : list
            charges and corresponding positions in angstroms as xyz coordinates

        Returns
        -------
        dict
            A dictionary with energy('energy') and gradient('gradients') information
             
        """

        topology, positions = self.convert_trajectory(traj)

        if charges is not None:
            self.set_external_charges(charges)

        info = self.compute_info(topology, positions, include_coulomb=include_coulomb, link_atoms=link_atoms, minimize=minimize)

        return info


    @abstractmethod
    def compute_info(self):
        pass

    @abstractmethod
    def set_external_charges(self):
        pass

    @abstractmethod
    def initialize(self):
        pass

    @abstractmethod
    def take_step(self, force):
        pass

    @abstractmethod
    def get_main_info(self):
        pass

    @abstractmethod
    def get_main_charges(self):
        pass

    @abstractmethod
    def convert_trajectory(self):
        pass

    @abstractmethod
    def convert_input(self):
        pass
    
    @abstractmethod
    def set_up_reporters(self):
        pass

    def get_qm_geometry(self):
        raise Exception('method not implemented for class')

    def build_qm_param(self):
        raise Exception('method not implemented for class')

    def optimize_geometry(self):
        raise Exception('method not implemented for class')
