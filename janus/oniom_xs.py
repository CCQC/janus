from .aqmmm import AQMMM
from .system import System
import numpy as np

class ONIOM_XS(AQMMM):

    def __init__(self, config, qm_wrapper, mm_wrapper):
        
        super().__init__(config, qm_wrapper, mm_wrapper, 'ONIOM-XS')

    def partition(self, qm_center=None, info=None): 
    
        if qm_center is None:
            qm_center = self.qm_center

        self.define_buffer_zone(qm_center)

        qm = System(qm_indices=self.qm_atoms, run_ID=self.run_ID, partition_ID='qm')

        self.systems[self.run_ID] = {}
        self.systems[self.run_ID][qm.partition_ID] = qm

        # the following only runs if there are groups in the buffer zone
        if self.buffer_groups:
            qm_bz = Partition(indices=self.qm_atoms, partition_ID='qm_bz')
            for key, value in self.buffer_groups.items():
                for idx in value:
                    qm_bz.qm_atoms.append(idx)

            # each partition has a copy of its buffer groups - 
            qm_bz.buffer_groups = self.buffer_groups

            self.systems[self.run_ID][qm_bz.partition_ID] = qm_bz

    def run_aqmmm(self):
        
        qm = self.systems[self.run_ID]['qm']

        if not self.buffer_groups:
            self.systems[self.run_ID]['qmmm_forces'] = qm.qmmm_energy
            self.systems[self.run_ID]['qmmm_energy'] = qm.qmmm_forces

        else:
            qm_bz = self.systems[self.run_ID]['qm_bz']
            lamda, d_lamda = self.get_switching_function(qm_bz)

            self.systems[self.run_ID]['qmmm_energy'] = \
            lamda*qm.qmmm_energy + (1-lamda)*qm_bz.qmmm_energy

            # needs work!
            print('need to add in d_lamda term for forces')
            forces = {}
            for f, coord in qm_bz.qmmm_forces.items():
                if f in qm.qmmm_forces:
                    forces[f] = (1-lamda)*coord + lamda*qm.qmmm_forces[f] 
                else: 
                    forces[f] = (1-lamda)*coord

            self.systems[self.run_ID]['qmmm_forces'] = forces


    def get_switching_function(self, partition):

        s = 0.0
        d_s = 0.0
        partition.switching_functions = []
        if partition.buffer_groups:
            for key, value in partition.buffer_groups.items():
                positions = self.get_qm_positions(value, as_string=False)
                COM = partition.compute_COM(positions)

                r_i = np.linalg.norm(COM - self.qm_center_xyz)
                s_i, d_s_i = self.compute_lamda_i(r_i)
                partition.switching_functions.append(s_i)
                s += s_i
                #d_s += d_s_i
                
            #print(" s",s)
            s *= 1/len(partition.switching_functions)
            print(partition.switching_functions)
        #d_s *= 1/len(partition.switching_functions)
        return s, d_s
            

    def compute_lamda_i(self, r_i):

        x_i = float((r_i - self.Rmin) / (self.Rmax - self.Rmin))

        
        lamda_i = 6*(x_i - 1/2)**5 - 5*(x_i - 1/2)**3 + (15/8)*(x_i - 1/2) + 1/2
        
        d_lamda_i = 30*(x_i - 1/2)**4 - 15*(x_i - 1/2)**2 + 15/8

        return lamda_i, d_lamda_i
        
