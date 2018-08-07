
class ONIOM_XS(AQMMM):

    def __init__(self):
        
        super().__init__(partition_scheme, trajectory=None)

    def partition(self, info): 

        self.define_buffer_zone()
        qm = Partition(indices=self.rmin_atoms, ID='qm')
        qm_bz = Partition(indices=self.rmin_atoms, ID='qm_bz')
        for key, value in self.buffer_groups.items():
            for idx in value:
                qm_bz.qm_atoms.append(idx)
    
        qm.positions = self.get_qm_positions(qm.qm_atoms)
        qm_bz.positions = self.get_qm_positions(qm_bz.qm_atoms)

        self.partitions[qm.ID] = qm 
        self.partitions[qm_bz.ID] = qm_bz

        return self.partitions

    def get_info(self):
        
        s = self.get_switching_function()
        qm = self.partitions['qm'] 
        qm_bz = self.partitions['qm_bz'] 
        self.energy = s*qm.energy + (1-s)*qm_bz.energy

    def get_switching_function(self):
        pass
        

