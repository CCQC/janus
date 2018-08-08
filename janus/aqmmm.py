from abc import ABC, abstractmethod
from .system import Partition
import mdtraj as md
"""
AQMMM class for adaptive QMMM computations
"""
class AQMMM(ABC):

    nm_to_angstrom = 10.0000000

    def __init__(self, partition_scheme, trajectory=None):

        self.partition_scheme = partition_scheme
        if self.partition_scheme == 'distance':
            # these values taken from original ONIOM-XS paper
            self.Rmin = 0.38 # in nm
            self.Rmax = 0.4 # in nm

        if trajectory is None:
            self.traj = md.load(system.mm_pdb_file)
        else:
            self.traj = trajectory

        # for now, need to define later
        self.qm_center = None
        # this needs to be np.array
        self.qm_center_xyz = None
        self.partitions = {}

    def save(self, ID, qmmm_forces, qmmm_energy):
        # find the appropriate partition object to save to
        self.partitions[ID].forces = qmmm_forces
        self.partitions[ID].energy = qmmm_energy


    def define_buffer_zone(self):

        if self.partition_scheme == 'distance': 
#            self.traj.xyz = positions
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


    def get_qm_positions(self, qm_atoms, as_string=True):
        """
        TODO:
        1. need to phase out getting qm_positions in the openmm wrapper
        2. need to phase out getting link atom stuff through the openmm wrapper
        - MDtraj can do ALL - just need to convert to openmm trajectory
        - this way would be more general and robust - the only thing is to make sure the qmmm 
         only things still work - not just with aqmmm
        """

        if as_string is True:
            out = ""
        if as_string is True:
            out = []
        line = '{:3} {: > 7.3f} {: > 7.3f} {: > 7.3f} \n '

        for idx in qm_atoms:
            x, y, z =   self.traj.xyz[0][idx][0]*AQMMM.nm_to_angstrom,\
                        self.traj.xyz[0][idx][1]*AQMMM.nm_to_angstrom,\
                        self.traj.xyz[0][idx][2]*AQMMM.nm_to_angstrom

            symbol = self.traj.topology.atom(idx).element.symbol
            
            if as_string is True:
                out += line.format(symbol, x, y, z)
            else:
                out.append([symbol, [x, y, z]])
            

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


    @abstractmethod
    def partition(self, info):
        pass

    @abstractmethod
    def get_info(self):
        pass
    







