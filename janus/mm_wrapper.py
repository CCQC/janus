"""
MM wrapper super class
"""
from abc import ABC, abstractmethod

class MM_wrapper(ABC):

    kjmol_to_au = 1/2625.5002 
    nm_to_angstrom = 10.0000000
    nm_to_bohr = 0.052917999999 
    kjmol_nm_to_au_bohr = kjmol_to_au*nm_to_bohr 
    au_bohr_to_kjmol_nm = 1/kjmol_nm_to_au_bohr

    def __init__(self, config, program):

        self.pdb_file = config['mm_pdb_file']

        if 'mm_temp' in config:
            self.temp = config['mm_temp']
        else:
            self.temp = float(300)

        if 'mm_step_size' in config:
            self.step_size = config['mm_step_size']
        else: 
            self.step_size = float(0.002)

        if 'embedding_method' in config:
            self.embedding_method = config['embedding_method']
        else:
            self.embedding_method = 'Mechanical'
        
        self.program = program


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

    
# still need to implement make zero function?
