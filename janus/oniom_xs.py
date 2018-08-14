from .aqmmm import AQMMM
from .system import Partition
import numpy as np

class ONIOM_XS(AQMMM):

    def __init__(self, partition_scheme, trajectory):
        
        super().__init__(partition_scheme, trajectory, 'ONIOM-XS')

    def partition(self, qm_center=None, info=None): 
    
        # need to first update everything with info!
        if info is not None:
            self.update_traj(info['positions'], info['topology'])

        if qm_center is None:
            qm_center = self.qm_center

        self.define_buffer_zone(qm_center)

        qm = Partition(indices=self.qm_atoms, ID='qm')
        qm.qm_positions = self.get_qm_positions(qm.qm_atoms)
        self.partitions[qm.ID] = qm 

        # the following only runs if there are groups in the buffer zone
        if self.buffer_groups:
            qm_bz = Partition(indices=self.qm_atoms, ID='qm_bz')
            for key, value in self.buffer_groups.items():
                for idx in value:
                    qm_bz.qm_atoms.append(idx)

                
            qm_bz.qm_positions = self.get_qm_positions(qm_bz.qm_atoms)
            # each partition has a copy of its buffer groups - 
            # good for later when there are multiple partitions with all different
            # buffer groups
            qm_bz.buffer_groups = self.buffer_groups

            self.partitions[qm_bz.ID] = qm_bz

        return self.partitions

    def get_info(self):
        
        qm = self.partitions['qm'] 

        if not self.buffer_groups:
            self.energy = qm.energy
            self.forces = qm.forces

        else:
            qm_bz = self.partitions['qm_bz'] 
            lamda, d_lamda = self.get_switching_function(qm_bz)
            self.energy = lamda*qm.energy + (1-lamda)*qm_bz.energy

            print('qm forces', qm.forces)
            print('qm_bz forces', qm_bz.forces)
            # needs work!
            print('need to add in d_lamda term for forces')
            self.forces = {}
            for f, coord in qm_bz.forces.items():
                if f in qm.forces:
                    self.forces[f] = (1-lamda)*coord + lamda*qm.forces[f] 
                else: 
                    self.forces[f] = (1-lamda)*coord
        
        

        return self.forces

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
        
        

