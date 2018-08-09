from .aqmmm import AQMMM
from copy import deepcopy

class HOT_SPOT(AQMMM):

    def __init__(self):
        
        super().__init__(partition_scheme, trajectory=None)

    def partition(self, info): 

        # need to first update everything with info!
    
        self.define_buffer_zone()

        qm = Partition(indices=self.rmin_atoms, ID='qm')
        qm.positions = self.get_qm_positions(qm.qm_atoms)
        self.partitions[qm.ID] = qm 

        # the following only runs if there are groups in the buffer zone
        if self.buffer_groups:
            qm_bz = Partition(indices=self.rmin_atoms, ID='qm_bz')
            for key, value in self.buffer_groups.items():
                for idx in value:
                    qm_bz.qm_atoms.append(idx)
        
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
            # does this create a deep copy??
            self.energy = qm.energy
            self.forces = qm.forces

        else:
            qm_bz = self.partitions['qm_bz'] 
            self.get_switching_function(qm_bz)
            
            # make sure this is deepcopied
            self.forces = deepcopy(qm_bz.forces)
            # counter for keeping track of lamda_i
            i = 0
            for key, value in self.buffer_groups.items():
                for idx in value: 
                    self.forces[idx] *= switching_functions[i]
                    self.forces[idx] += (1 - switching_functions[i])*qm.forces[idx]
                i += 1

    def get_switching_function(self, partition):

        for key, value in self.buffer_groups.items():
            positions = self.get_qm_positions(value, as_string=False)
            COM = partition.compute_COM(positions)

            r_i = COM - self.qm_center_xyz
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
