import mdtraj as md
"""
AQMMM class for adaptive QMMM computations
"""
class AQMMM(object):

    nm_to_angstrom = 10.0000000
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

        # for adding identifying water buffer groups
        groups = {}
        top = self.traj.topology

        for i in self.buffer_atoms:
            if top.atom(i).residue.is_water:
                idx = top.atom(i).residue.index
                if idx not in groups.keys():
                    groups[idx] = []
                    for a in top.residue(idx).atoms:
                        if a.index not in groups[idx]:
                            group[idx].append(a.index)

        self.buffer_groups = groups


    def oniom_xs(partition=False, get_info=False):
        
        if partition is True:
            qm = Partition(indices=self.rmin_atoms, ID=1)
            qm_bz = Partition(indices=self.rmin_atoms, ID=1)
            for key, value in self.buffer_groups.items():
                for idx in value:
                    qm_bz.indicies.append(idx)
        
            qm.positions = self.get_qm_positions(qm.indicies)
            qm_bz.positions = self.get_qm_positions(qm_bz.indicies)

            self.partitions = [qm, qm_bz]

            return self.partitions

        if get_info is True:
            pass

    def get_qm_positions(self, qm_atoms):
        """
        TODO:
        1. need to phase out getting qm_positions in the openmm wrapper
        2. need to phase out getting link atom stuff through the openmm wrapper
        - MDtraj can do ALL - just need to convert to openmm trajectory
        - this way would be more general and robust - the only thing is to make sure the qmmm 
         only things still work - not just with aqmmm
        """
        out = ""
        line = '{:3} {: > 7.3f} {: > 7.3f} {: > 7.3f} \n '

        for idx in qm_atoms:
            for atom in self.traj.topology.atoms():
                if atom.index == idx:
                    x, y, z =   self.traj.xyz[0][idx][0]*AQMMM.nm_to_angstrom,\
                                self.traj.xyz[0][idx][1]*AQMMM.nm_to_angstrom,\
                                self.traj.xyz[0][idx][2]*AQMMM.nm_to_angstrom
                    out += line.format(atom.element.symbol, x, y, z)

        ## if there are bonds that need to be cut
        #if self._boundary_bonds:
        #    # Need to add if statement for any treatments that don't need link atoms
        #    #if self._system.boundary_treatment !=
        #    for atom in self.link_atoms:
        #        pos = self.link_atoms[atom]['link_positions']*MM_wrapper.nm_to_angstrom
        #        x, y, z = pos[0], pos[1], pos[2]
        #        out += line.format(self.link_atoms[atom]['link_atom'], x, y, z)
        return out



    def get_Rmin(self):
        return self.Rmin

    def get_Rmax(self):
        return self.Rmax

    def set_Rmin(self, Rmin):
        self.Rmin = Rmin

    def set_Rmax(self, Rmax):
        self.Rmax = Rmax
    


class Partition(object):
    def __init__(self, indices, ID):
        self.indices = indices
        self.ID = ID
        self.qm_positions = None






