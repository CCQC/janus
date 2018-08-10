from abc import ABC, abstractmethod
from copy import deepcopy
import mdtraj as md
import numpy as np

"""
AQMMM class for adaptive QMMM computations
"""
class AQMMM(ABC):

    nm_to_angstrom = 10.0000000

    def __init__(self, partition_scheme, trajectory, scheme=None):

        self.partition_scheme = partition_scheme
        if self.partition_scheme == 'distance':
            # from oniom-xs paper values are 0.38 and 0.4
            self.Rmin = 0.38 # in nm
            self.Rmax = 0.45 # in nm

        self.traj = md.load(trajectory)

        # for now, need to define later
        #self.qm_center = None
        # this needs to be np.array
        self.partitions = {}

    def save(self, ID, qmmm_forces, qmmm_energy):
        # find the appropriate partition object to save to
        self.partitions[ID].forces = qmmm_forces
        self.partitions[ID].energy = qmmm_energy


    def define_buffer_zone(self, qm_center):
        # qm_center needs to be in list form
        self.qm_center = qm_center
        self.qm_center_xyz = self.traj.xyz[0][qm_center[0]]

        if self.partition_scheme == 'distance': 
#            self.traj.xyz = positions
            rmin_atoms = md.compute_neighbors(self.traj, self.Rmin, qm_center)
            rmax_atoms = md.compute_neighbors(self.traj, self.Rmax, qm_center)
            self.buffer_atoms = np.setdiff1d(rmax_atoms, rmin_atoms)
            self.qm_atoms = rmin_atoms[0].tolist()
            self.qm_atoms.append(qm_center[0])
            self.qm_atoms.sort()
            

        # for adding identifying water buffer groups
        groups = {}
        top = self.traj.topology

        for i in self.buffer_atoms:
            # since if a hydrogen is in buffer zone with link atoms the qm would be the same as qm_bz, and the center 
             # of mass would not be in the buffer zone
            # only if oxygen in buffer zone
            if (top.atom(i).residue.is_water and top.atom(i).element.symbol == 'O'):
                idx = top.atom(i).residue.index
                if idx not in groups.keys():
                    groups[idx] = []
                    for a in top.residue(idx).atoms:
                        if a.index not in groups[idx]:
                            groups[idx].append(a.index)
                        # if oxygen in buffer zone remove any hydrogens in qm area from definition of qm atoms
                        if a.index in self.qm_atoms:
                            self.qm_atoms.remove(a.index)

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
        if as_string is False:
            out = []
        line = '{:3} {: > 7.3f} {: > 7.3f} {: > 7.3f} \n '

        for idx in qm_atoms:
            x, y, z =   self.traj.xyz[0][idx][0],\
                        self.traj.xyz[0][idx][1],\
                        self.traj.xyz[0][idx][2]

            symbol = self.traj.topology.atom(idx).element.symbol
            
            if as_string is True:
                # convert to angstroms
                out += line.format(symbol, x*10, y*10, z*10)

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
    







