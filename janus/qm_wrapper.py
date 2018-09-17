"""
QM wrapper super class
"""
from abc import ABC, abstractmethod
class QM_wrapper(ABC):

    def __init__(self, param, program):
        self.program = program
        self.param = param

        self.qm_param = None
        self.external_charges = None
        self.charges = None

    @abstractmethod
    def build_qm_param(self):
        pass

    def run_qm(self, geometry, total_elec):
        """
        Gets the energy and gradient from a QM computation of the primary subsystem 
        NEED TO RENAME

        Parameters
        ----------
        None

        Returns
        -------
        A dictionary with energy and gradient information
        """
        self.set_qm_geometry(geometry, total_elec)
        if not self.qm_param:
            self.build_qm_param()
        self.compute_energy_and_gradient()

        self.info = {}
        self.info['energy'] = self.energy
        self.info['gradients'] = self.gradient
        
        return self.info

    def set_qm_geometry(self, qm_geometry, total_elec):

        self.qm_geometry = qm_geometry
        self.total_elec = total_elec

        if total_elec % 2 != 0:
            self.total_elec += self.charge
            if total_elec % 2 != 0:
                self.is_close_shelled = False
            
    def set_external_charges(self, charges):
        
        self.external_charges = charges

