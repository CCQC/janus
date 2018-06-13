"""
QM wrapper super class
"""
from abc import ABC, abstractmethod
class QM_wrapper(ABC):

    def __init__(self, system, program):

        self._system = system
        self._program = program
        self._qm = {}
    
    @abstractmethod
    def qm_info(self):
        pass

    def get_qm(self):
        """
        Gets the energy and gradient from a QM computation of the primary subsystem 

        Parameters
        ----------
        None

        Returns
        -------
        A dictionary with energy and gradient information
        """
        if not self._qm:
            self.qm_info()
        return self._qm 
