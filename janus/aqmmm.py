import mdtraj as md
"""
AQMMM class for adaptive QMMM computations
"""
class AQMMM(object):

    def __init__(self, system, trajectory=None):

        self.scheme = system.aqmmm_scheme
        self.partition_scheme = system.aqmmm_partition_scheme
        if self.partition_scheme == 'distance':
            self.Rmin = 0.5 # in nm
            self.Rmax = 0.6 # in nm

        if trajectory is None:
            self.traj = md.load(system.mm_pdb_file)
        else:
            self.traj = trajectory

        # for now, need to define later
        self.qm_center = None

    def partition(self, info):

        self.define_buffer_zone(position)

        # make this into class structure?
        if self.scheme == 'ONIOM-XS':
            self.oniom_xs(partition=True)

    def save(self):
        pass

    def get_info(self):

        if self.scheme == 'ONIOM-XS':
            self.oniom_xs(get_info=True)
        
        return self.info

    def define_buffer_zone(self, positions):

        if self.partition_scheme == 'distance': 
            self.traj.xyz = positions
            self.rmin_atoms = md.compute_neighbors(self.traj, self.Rmin, self.qm_center)
            self.rmax_atoms = md.compute_neighbors(self.traj, self.Rmax, self.qm_center)
            self.buffer_atoms = np.setdiff1d(self.rmax_atoms, self.rmin_atoms)
            

    def oniom_xs(partition=False, get_info=False):
        
        if partition is True:
            pass

        if get_info is True:
            pass


    def get_Rmin(self):
        return self.Rmin

    def get_Rmax(self):
        return self.Rmax

    def set_Rmin(self, Rmin):
        self.Rmin = Rmin

    def set_Rmax(self, Rmax):
        self.Rmax = Rmax

