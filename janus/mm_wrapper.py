from abc import ABC, abstractmethod

class MM_wrapper(ABC):
    """
    A super class for all molecular mechanics wrappers
    """

    kjmol_to_au = 1/2625.5002 
    nm_to_angstrom = 10.0000000
    nm_to_bohr = 0.052917999999 
    kjmol_nm_to_au_bohr = kjmol_to_au*nm_to_bohr 
    au_bohr_to_kjmol_nm = 1/kjmol_nm_to_au_bohr

    def __init__(self, param, program):

        if param['mm_pdb_file']
            self.pdb_file = param['mm_pdb_file']
        self.param = param
        self.program = program
        self.main_info = None
        self.main_charges = None

        super().__init__()

    def get_energy_and_gradient(self, traj, include_coulomb='all', link_atoms=None, minimize=False):

        topology, positions = self.convert_trajectory(traj)

        if charges is not None:
            self.set_external_charges(charges)

        info = self.mm_wrapper.compute_mm(topology, positions, include_coulomb=include_coulomb, link_atoms=link_atoms, minimize=minimize)

        return info

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
    def compute_mm(self):
        pass

    @abstractmethod
    def get_main_charges(self):
        pass

    @abstract_method
    def convert_trajectory(self):
        pass

    @abstract_method
    def equilibrate(self):
        pass
    
