from .aqmmm import AQMMM
from .system import System
import numpy as np
from copy import deepcopy

class HotSpot(AQMMM):

    def __init__(self, config, qm_wrapper, mm_wrapper):
        
        super().__init__(config, qm_wrapper, mm_wrapper, 'Hot-Spot')

    def partition(self, qm_center=None, info=None): 
    
        if qm_center is None:
            qm_center = self.qm_center

        self.define_buffer_zone(qm_center)

        qm = System(qm_indices=self.qm_atoms, run_ID=self.run_ID, partition_ID='qm')

        # the following only runs if there are groups in the buffer zone
        if self.buffer_groups:

            qm.buffer_groups = deepcopy(self.buffer_groups)

            for key, value in self.buffer_groups.items():
                for idx in value:
                    qm.qm_atoms.append(idx)
                
        self.systems[self.run_ID] = {}
        self.systems[self.run_ID][qm.partition_ID] = qm

    def run_aqmmm(self):
        
        qm = self.systems[self.run_ID]['qm']

        # do I need to do deepcopy?
        if not self.buffer_groups:
            self.systems[self.run_ID]['qmmm_forces'] = qm.qmmm_energy
            self.systems[self.run_ID]['qmmm_energy'] = qm.qmmm_forces

        else:

            forces = deepcopy(qm.forces)
            for key, value in self.buffer_groups.items():
                for idx in value: 
                    forces[idx] *= self.buffer_switching_functions[key][0]

            self.systems[self.run_ID]['qmmm_forces'] = forces

    def compute_lamda_i(self, r_i):

        if r_i <= self.Rmin:
            lamda_i = 1


        elif r_i > self.Rmax:
            lamda_i = 0

        else:
            lamda_i = (self.Rmax**2 - r_i**2)**2 * (self.Rmax**2 + 2*(r_i**2) - 3*(self.Rmin**2))
            lamda_i *= 1/(self.Rmax**2 - self.Rmin**2)**3

        return lamda_i, None
