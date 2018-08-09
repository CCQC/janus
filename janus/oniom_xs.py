from .aqmmm import AQMMM

class ONIOM_XS(AQMMM):

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
            self.energy = qm.energy
            self.forces = qm.forces

        else:
            qm_bz = self.partitions['qm_bz'] 
            lamda, d_lamda = self.get_switching_function(qm_bz)
            self.energy = lamda*qm.energy + (1-lamda)*qm_bz.energy
            # needs work!
            #self.forces = lamda*qm.forces + (1-lamda)*qm.forces + d_lamda*(qm.energy - qm_bz.energy)

        return self.forces

    def get_switching_function(self, partition):

        s = 0.0
        #d_s = 0.0
        for key, value in self.buffer_groups.items():
            positions = self.get_qm_positions(value, as_string=False)
            COM = partition.compute_COM(positions)

            r_i = COM - self.qm_center_xyz
            s_i, d_s_i = self.compute_s_i(r_i)
            partition.switching_functions.append(s_i)
            s += s_i
            #d_s += d_s_i
            
        s *= 1/len(partition.switching_functions)
        #d_s *= 1/len(partition.switching_functions)
        return s #, d_s
            

    def compute_s_i(self, r_i):

        x_i = (r_i - self.Rmin) / (self.Rmax - self.Rmin)
        
        lamda_i = 6*(x_i - 1/2)**5 - 5(x_i - 1/2)**3 + (15/8)*(x_i - 1/2) + 1/2
        
        d_lamda_i = 30*(x_i - 1/2)**4 - 15(x_i - 1/2)**2 + 15/8

        return s_i, d_s_i
        
        

