from abc import ABC, abstractmethod
from copy import deepcopy
import mdtraj as md
import numpy as np
import mendeleev as mdlv
from .qmmm import QMMM
from .system import Buffer

"""
AQMMM class for adaptive QMMM computations
"""
class AQMMM(ABC, QMMM):

    nm_to_angstrom = 10.0000000

    def __init__(self, config, qm_wrapper, mm_wrapper):
        super().__init__(config, qm_wrapper, mm_wrapper)

        self.aqmmm_scheme = config['aqmmm_scheme']
        self.partition_scheme = config['aqmmm_partition_scheme']
        self.Rmin = config['Rmin']
        self.Rmax = config['Rmax']
        # do not include options of computing the qm center with the program - 
        # might need this functionality later
        self.qm_center = config['qm_center']

        self.buffer_groups = {}

    def run_qmmm(self,main_info):

        self.update_traj(main_info['positions'], main_info['topology'])
        self.partition()
            
        for i, system in self.systems[self.run_ID].items():

            self.qm_atoms = deepcopy(system.qm_atoms)

            if self.embedding_method =='Mechanical':
                self.mechanical(system, main_info)
            elif self.embedding_method =='Electrostatic':
                self.electrostatic(system, main_info)
            else:
                print('only mechanical and electrostatic embedding schemes implemented at this time')

        self.run_aqmmm()
        self.run_ID += 1

    def define_buffer_zone(self, qm_center):
        # qm_center needs to be in list form
        self.qm_center = qm_center
        self.qm_center_xyz = self.traj.xyz[0][qm_center]

        if self.partition_scheme == 'distance': 
#            self.traj.xyz = positions
            rmin_atoms = md.compute_neighbors(self.traj, self.Rmin, qm_center)
            rmax_atoms = md.compute_neighbors(self.traj, self.Rmax, qm_center)
            self.buffer_atoms = np.setdiff1d(rmax_atoms, rmin_atoms)
            self.qm_atoms = rmin_atoms[0].tolist()
            self.qm_atoms.append(qm_center[0])
        
        
        self.edit_qm_atoms()

        # for adding identifying water buffer groups
        top = self.topology

        for i in self.buffer_atoms:
            # since if a hydrogen is in buffer zone with link atoms the qm would be the same as qm_bz, and the center 
             # of mass would not be in the buffer zone
            # only if oxygen in buffer zone
            if (top.atom(i).residue.is_water and top.atom(i).element.symbol == 'O'):
                idx = top.atom(i).residue.index
                if idx not in self.buffer_groups.keys():
                    buf = Buffer(ID=idx)
                    for a in top.residue(idx).atoms:
                        if a.index not in buf.atoms:
                            buf.atoms.append(a.index)

                    self.buffer_groups[idx] = buf

        if self.buffer_groups:
            self.get_buffer_info()

    def get_buffer_info(self):

        self.buffer_distance = {}

        for i, buf in self.buffer_groups.items():

            xyz, buf.atom_weights, buf.weight_ratio = self.compute_COM(buf.atoms)
            buf.COM_coord = np.array(xyz)
            buf.r_i = np.linalg.norm(buf.COM_coord - self.qm_center_xyz)
            self.buffer_distance[i] = buf.r_i
            buf.s_i, buf.d_s_i = self.compute_lamda_i(buf.r_i)

    def compute_COM(self, atoms):
        
        xyz = np.zeros(3)
        M = 0

        atom_weight = {}
        weight_ratio = {}
        for i in atoms:

            symbol = self.traj.topology.atom(i).element.symbol
            m = mdlv.element(symbol).atomic_weight
            # this gives positions in nm
            position =  np.array(self.traj.xyz[0][i])

            M += m
            xyz += m * position
            atom_weight[i] = m

        for i in atoms:
            weight_ratio[i] = atom_weight[i]/M
            
        xyz *= 1/M
        
        return xyz, atom_weight, weight_ratio

    def compute_lamda_i(self, r_i):

        x_i = float((r_i - self.Rmin) / (self.Rmax - self.Rmin))

        lamda_i = -6*((x_i)**5) + 15*((x_i)**4) - 10*((x_i)**3) + 1

        d_lamda_i = (-30*(x_i)**4  + 60*(x_i)**3 - 30*(x_i)**2)
        d_lamda_i *= 1/(r_i * (self.Rmax - self.Rmin))

        return lamda_i, d_lamda_i

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
    def run_aqmmm(self):
        pass


