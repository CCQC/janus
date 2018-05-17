"""
MM wrapper super clasecond_subsys
"""
from abc import ABC, abstractmethod

class MM_wrapper(ABC):

    kjmol_to_au = 1/2625.5 
    nm_to_angstrom = 10.0

    def __init__(self, system, program):
        
        self._program = program
        self._system = system
        self._second_subsys = {}
        self._primary_subsys = {}
        self._entire_sys = {}
        self._boundary = {} 
        self._boundary['energy'] = None
        self._qm_positions = None

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

#    @abstractmethod
#    def make_zero_energy(self):
#        pass

    def get_second_subsys(self):
        """
        mm energy and gradients on second_subsys only
        """
        if not self._second_subsys:
            self.second_subsys_info()
        return self._second_subsys

    def get_primary_subsys(self, link=False):
        """
        mm energy and gradients on primary_subsys only
        """
        if not self._primary_subsys:
            self.primary_subsys_info(link)
        return self._primary_subsys
    
    def get_entire_sys(self):
        """
        mm energy and gradients on entire_sys 
        """
        if not self._entire_sys:
            self.entire_sys_info()
        return self._entire_sys

    def get_boundary(self):
        """
        energy and gradients for interaction between 
        """
        if self._boundary['energy'] is None:
            self.boundary_info()
        return self._boundary

    def get_qm_positions(self):
        if self._qm_positions is None:
            self.qm_positions()
        return self._qm_positions


    
