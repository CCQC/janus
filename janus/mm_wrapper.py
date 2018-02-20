"""
MM wrapper super class
"""
from abc import ABC, abstractmethod

class MM_wrapper(ABC):

    kjmol_to_au = 1/2625.5 
    nm_to_angstrom = 10.0

    def __init__(self, system, program):
        
        self._program = program
        self._system = system
        self._ss = None
        self._ps = None
        self._es = None
        self._ps_ss = None
        self._ps_ss = None
        self._qm_positions = None

        super().__init__()

    @abstractmethod
    def ss_info(self):
        pass

    @abstractmethod
    def ps_info(self):
        pass

    @abstractmethod
    def es_info(self):
        pass

    @abstractmethod
    def ps_ss_info(self):
        pass

    @abstractmethod
    def qm_positions(self):
        pass

    def get_ss(self):
        """
        mm energy and gradients on ss only
        """
        if self._ss is None:
            self.ss_info()
        return self._ss

    def get_ps(self):
        """
        mm energy and gradients on ps only
        """
        if self._ps is None:
            self.ps_info()
        return self._ps
    
    def get_es(self):
        """
        mm energy and gradients on es 
        """
        if self._es is None:
            self.es_info()
        return self._es

    def get_ps_ss(self):
        """
        energy and gradients for interaction between 
        """
        if self._ps_ss is None:
            self.ps_ss_info()
        return self._ps_ss

    def get_qm_positions(self):
        if self._qm_positions is None:
            self.qm_positions()
        return self._qm_positions

    def make_zero_energy(self):
        pass

    
