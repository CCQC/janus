from .aqmmm import AQMMM
from .system import Partition
import numpy as np

class ONIOM_XS(AQMMM):

    def __init__(self, partition_scheme, trajectory):
        
        super().__init__(partition_scheme, trajectory, 'ONIOM-XS')

    def partition(self, qm_center, info=None): 
    
        # need to first update everything with info!
        if info is not None:
            self.traj = info.traj

        self.define_buffer_zone(qm_center)

        qm = Partition(indices=self.qm_atoms, ID='qm')
        qm.positions = self.get_qm_positions(qm.qm_atoms)
        self.partitions[qm.ID] = qm 
        print(qm.qm_atoms)

        # the following only runs if there are groups in the buffer zone
        if self.buffer_groups:
            qm_bz = Partition(indices=self.qm_atoms, ID='qm_bz')
            for key, value in self.buffer_groups.items():
                for idx in value:
                    if idx not in qm_bz.qm_atoms:
                        qm_bz.qm_atoms.append(idx)

            print(qm.qm_atoms)
            print(qm_bz.qm_atoms)
                
            qm_bz.positions = self.get_qm_positions(qm_bz.qm_atoms)
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
            lamda = self.get_switching_function(qm_bz)
            print('qm_energy', qm.energy)
            print('qm_energy', qm_bz.energy)
            print('lamda', lamda)
            self.energy = lamda*qm.energy + (1-lamda)*qm_bz.energy
            # needs work!
            #self.forces = lamda*qm.forces + (1-lamda)*qm.forces + d_lamda*(qm.energy - qm_bz.energy)
            self.forces = None
        print(self.energy)

        return self.forces

    def get_switching_function(self, partition):

        s = 0.0
        #d_s = 0.0
        if partition.buffer_groups:
            for key, value in partition.buffer_groups.items():
                print('value', value)
                positions = self.get_qm_positions(value, as_string=False)
                COM = partition.compute_COM(positions)
                print(COM)
                print(self.qm_center_xyz)

                r_i = np.linalg.norm(COM - self.qm_center_xyz)
                print(r_i)
                s_i, d_s_i = self.compute_lamda_i(r_i)
                partition.switching_functions.append(s_i)
                s += s_i
                #d_s += d_s_i
                
                print("current s",s)
            print(" s",s)
            s *= 1/len(partition.switching_functions)
        #d_s *= 1/len(partition.switching_functions)
        return s #, d_s
            

    def compute_lamda_i(self, r_i):

        x_i = float((r_i - self.Rmin) / (self.Rmax - self.Rmin))

        
        lamda_i = 6*(x_i - 1/2)**5 - 5*(x_i - 1/2)**3 + (15/8)*(x_i - 1/2) + 1/2
        
        d_lamda_i = 30*(x_i - 1/2)**4 - 15*(x_i - 1/2)**2 + 15/8

        return lamda_i, d_lamda_i
        
        

