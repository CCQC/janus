from .aqmmm import AQMMM
from .system import Partition
import numpy as np
from copy import deepcopy

class HotSpot(AQMMM):

    def __init__(self, partition_scheme, trajectory):
        
        super().__init__(partition_scheme, trajectory, 'Hot-Spot')

    def partition(self, qm_center=None, info=None): 
    
        # need to first update everything with info!
        if info is not None:
            self.update_traj(info['positions'], info['topology'])

        if qm_center is None:
            qm_center = self.qm_center

        self.define_buffer_zone(qm_center)

        qm = Partition(indices=self.qm_atoms, ID='qm')

        # the following only runs if there are groups in the buffer zone
        if self.buffer_groups:

            qm.buffer_groups = self.buffer_groups

            for key, value in self.buffer_groups.items():
                for idx in value:
                    qm.qm_atoms.append(idx)
                
        qm.positions = self.get_qm_positions(qm.qm_atoms)
        self.partitions[qm.ID] = qm 

        return self.partitions

    def get_info(self):
        
        qm = self.partitions['qm'] 

        if not self.buffer_groups:
            # does this create a deep copy??
            self.energy = qm.energy
            self.forces = qm.forces

        else:
            self.get_switching_function(qm)
            
            # make sure this is deepcopied
            self.forces = deepcopy(qm.forces)
            print(self.forces)
            # counter for keeping track of lamda_i
            i = 0
            for key, value in self.buffer_groups.items():
                for idx in value: 
                    self.forces[idx] *= qm.switching_functions[i]
                i += 1

        return self.forces

    def get_switching_function(self, partition):

        partition.switching_functions = []
        if partition.buffer_groups:
            for key, value in self.buffer_groups.items():
                positions = self.get_qm_positions(value, as_string=False)
                COM = partition.compute_COM(positions)

                r_i = np.linalg.norm(COM - self.qm_center_xyz)
                lamda_i = self.compute_lamda_i(r_i)
                partition.switching_functions.append(lamda_i)
            

    def compute_lamda_i(self, r_i):

        if r_i <= self.Rmin:
            lamda_i = 1

        elif r_i > self.Rmax:
            lamda_i = 0

        else:
            lamda_i = (self.Rmax**2 - r_i**2)**2 * (self.Rmax**2 + 2*r_i**2 - 3*self.Rmin**2)
            lamda_i *= 1/(self.Rmax**2 - self.Rmin**2)**2

        return lamda_i
