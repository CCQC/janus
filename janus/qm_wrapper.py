from abc import ABC, abstractmethod
"""
QM wrapper super class
"""

class QM_wrapper(ABC):

    def __init__(self, param, program):
        self.program = program
        self.param = param

        self.qm_param = None
        self.external_charges = None
        self.charges = None
        self.is_open_shelled = False

    @abstractmethod
    def build_qm_param(self):
        pass

    def run_qm(self, geometry, total_elec, minimize=False):
        """
        Gets the energy and gradient from a QM computation of the primary subsystem 

        Parameters
        ----------
        geometry: A string with the geometry information 
                  for the QM region
        total_elec: An int of the total number of electrons 
                    present in the QM region

        Returns
        -------
        A dictionary with energy('energy') and gradient('gradients') information

        Examples
        --------
        run_qm(geom, 10)
        """

        self.set_qm_geometry(geometry, total_elec)
        if not self.qm_param:
            self.build_qm_param()

        if minimize is True:
            geom = self.optimize_geometry()
        else:
            self.compute_energy_and_gradient()

        self.info = {}
        self.info['energy'] = self.energy
        self.info['gradients'] = self.gradient
        
        return self.info

    def set_qm_geometry(self, qm_geometry, total_elec):
        """
        Sets the geometry for the QM region as 
        self.qm_geometry. Also determines if the system is open shelled
        based on the number of electrons

        Parameters
        ----------
        geometry: A string with the geometry information 
                  for the QM region
        total_elec: An int of the total number of electrons 
                    present in the QM region

        Returns
        -------
        None
        
        Examples
        --------
        set_qm_geometry(geom, 10)
        """

        self.qm_geometry = qm_geometry
        self.total_elec = total_elec

        if self.total_elec % 2 != 0:
            self.total_elec += self.charge   # takes charge into account
            if self.total_elec % 2 != 0:
                self.is_open_shelled = True
            
    def set_external_charges(self, charges):
        """
        Sets the charges due to the MM point charges
        as self.external_charges

        Parameters
        ----------
        charges: a list with the charge and position(in angstroms) for all 
                 particles outside the QM region 

        Returns
        -------
        None

        Examples
        --------
        set_external_charges([[-.45, 0.0, 0.0, 0.0]])
        set_external_charges(charges)
        """
        
        self.external_charges = charges

