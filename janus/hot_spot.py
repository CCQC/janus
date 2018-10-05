from .aqmmm import AQMMM
from .system import System
import numpy as np
from copy import deepcopy

class HotSpot(AQMMM):
    """
    Class for the Hot-Spot adaptive QM/MM method.
    Inherits from AQMMM class
    """

    def __init__(self, param, hl_wrapper, ll_wrapper, md_simulation_program):
        """
        Initializes the HotSpot class object
    
        Parameters
        ----------
        See parameters for AQMMM class 

        """
        
        super().__init__(param, hl_wrapper, ll_wrapper, md_simulation_program, 'Hot-Spot')

    def partition(self, qm_center=None): 
        """
        Finds the partitions as required by the Hot-Spot method 
        and saves each partition as a system object.
        Saves all systems in the dictionary self.systems

        Parameters
        ----------
        qm_center : list 
            atoms that define the qm center, default is None

        """
    
        if qm_center is None:
            qm_center = self.qm_center

        self.define_buffer_zone(qm_center)

        qm = System(qm_indices=self.qm_atoms, qm_residues=self.qm_residues, run_ID=self.run_ID, partition_ID='qm')

        # the following only runs if there are groups in the buffer zone
        if self.buffer_groups:

            qm.buffer_groups = deepcopy(self.buffer_groups)

            for i, buf in self.buffer_groups.items():
                qm.qm_residues.append(i)
                for idx in buf.atoms:
                    qm.qm_atoms.append(idx)
                
        self.systems[self.run_ID] = {}
        self.systems[self.run_ID][qm.partition_ID] = qm

    def run_aqmmm(self):
        """
        Interpolates the energy and gradients from each partition
        according to the Hot-Spot method
        """
        
        qm = self.systems[self.run_ID]['qm']

        # do I need to do deepcopy?
        if not self.buffer_groups:
            self.systems[self.run_ID]['qmmm_energy'] = qm.qmmm_energy
            self.systems[self.run_ID]['qmmm_forces'] = qm.qmmm_forces

        else:

            forces = deepcopy(qm.qmmm_forces)
            for i, buf in self.buffer_groups.items():
                for idx in buf.atoms:
                    forces[idx] *= buf.s_i

            self.systems[self.run_ID]['qmmm_forces'] = forces

    def compute_lamda_i(self, r_i):
        """
        Computes the switching function of the Hot-Spot method
        and overrides the compute_lamda_i function from the AQMMM class
        
        Parameters
        ----------
        r_i : float 
            the distance between the qm center and the COM in angstroms 

        Returns
        -------
        float
            lamda_i, unitless
        None

        """

        if r_i <= self.Rmin:
            lamda_i = 1


        elif r_i > self.Rmax:
            lamda_i = 0

        else:
            lamda_i = (self.Rmax**2 - r_i**2)**2 * (self.Rmax**2 + 2*(r_i**2) - 3*(self.Rmin**2))
            lamda_i *= 1/(self.Rmax**2 - self.Rmin**2)**3

        return lamda_i, None
