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

        self.pdb_file = param['mm_pdb_file']
        self.param = param
        self.program = program
        self.main_info = None
        self.main_charges = None

        super().__init__()


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

    
    
# still need to implement make zero function?
