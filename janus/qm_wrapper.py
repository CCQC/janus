"""
QM wrapper super class
"""
from abc import ABC, abstractmethod
class QM_wrapper(ABC):

    def __init__(self, system, program):

        self._system = system
        self._program = program
        self._qm = None
    
    @abstractmethod
    def qm_info(self):
        print('using qm wrapper')

    def get_qm(self):
        if self._qm is None:
            self.qm_info()
        return self._qm 
