from .aqmmm import AQMMM
from .system import System
import numpy as np

class ONIOM_XS(AQMMM):

    def __init__(self, config, qm_wrapper, mm_wrapper):
        
        super().__init__(config, qm_wrapper, mm_wrapper)

    def partition(self, qm_center=None, info=None): 
    
        if qm_center is None:
            qm_center = self.qm_center

        self.define_buffer_zone(qm_center)

        qm = System(qm_indices=self.qm_atoms, run_ID=self.run_ID, partition_ID='qm')

        self.systems[self.run_ID] = {}
        self.systems[self.run_ID][qm.partition_ID] = qm

        # the following only runs if there are groups in the buffer zone
        if self.buffer_groups:
            qm_bz = System(qm_indices=self.qm_atoms, run_ID=self.run_ID, partition_ID='qm_bz')
            for i, buf in self.buffer_groups.items():
                for idx in buf.atoms:
                    qm_bz.qm_atoms.append(idx)

            # each partition has a copy of its buffer groups - 
            qm_bz.buffer_groups = self.buffer_groups

            self.systems[self.run_ID][qm_bz.partition_ID] = qm_bz

    def run_aqmmm(self):
        
        qm = self.systems[self.run_ID]['qm']

        if not self.buffer_groups:
            self.systems[self.run_ID]['qmmm_energy'] = qm.qmmm_energy
            self.systems[self.run_ID]['qmmm_forces'] = qm.qmmm_forces

        else:
            qm_bz = self.systems[self.run_ID]['qm_bz']
            lamda = self.get_switching_function()
            

            self.systems[self.run_ID]['qmmm_energy'] = \
            (1- lamda)*qm.qmmm_energy + lamda*qm_bz.qmmm_energy

            # needs work!
            print('need to add in d_lamda term for forces')
            forces = {}
            for f, coord in qm_bz.qmmm_forces.items():
                if f in qm.qmmm_forces:
                    forces[f] = lamda*coord + (1-lamda)*qm.qmmm_forces[f] 
                else: 
                    forces[f] = lamda*coord

            self.systems[self.run_ID]['qmmm_forces'] = forces


    def get_switching_function(self):

        s = 0.0
                
        for i, buf in self.buffer_groups.items():
            s += buf.s_i
            
        s *= 1/len(self.buffer_groups)

        return s
            

        
