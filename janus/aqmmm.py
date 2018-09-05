from abc import ABC, abstractmethod
from copy import deepcopy
import mdtraj as md
import numpy as np
import mendeleev as mdlv
from .qmmm import QMMM

"""
AQMMM class for adaptive QMMM computations
"""
class AQMMM(ABC, QMMM):

    nm_to_angstrom = 10.0000000

    def __init__(self, config, qm_wrapper, mm_wrapper):
        super().__init__(config, qm_wrapper, mm_wrapper)


        # for now, need to define later
        #self.qm_center = None
        # this needs to be np.array

        if 'aqmmm_scheme' in config:
            self.aqmmm_scheme = config['aqmmm_scheme']
        else: 
            self.aqmmm_scheme = 'ONIOM-XS'

        if 'aqmmm_partition_scheme' in config:
            self.partition_scheme = config['aqmmm_partition_scheme']
        else:
            self.partition_scheme = 'distance'

        if (self.partition_scheme == 'distance' and 'Rmin' in config):
            # from oniom-xs paper 0.38
            self.Rmin = config['Rmin']
        else:
            self.Rmin = 0.38 # in nm

        if (self.partition_scheme == 'distance' and 'Rmax' in config):
            # from oniom-xs paper 0.4
            self.Rmax = config['Rmax']
        else:
            self.Rmax = 0.45 
        
        if 'qm_center' in config:
            self.qm_center = config['qm_center']
        else:
        # this does not include options of computing the qm center with the program - 
        # might need this functionality later
            self.qm_center = [0]

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
        groups = {}
        top = self.topology

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

        self.buffer_groups = groups
        if groups:
            self.get_buffer_info()

    def get_buffer_info(self):

        self.buffer_switching_functions = {}
        self.buffer_distance = {}

        for key, value in self.buffer_groups.items():

            COM = self.compute_COM(value)
            r_i = np.linalg.norm(COM - self.qm_center_xyz)
            self.buffer_distance[key] = r_i
            s_i, d_s_i = self.compute_lamda_i(r_i)
            self.buffer_switching_functions[key] = [s_i, d_s_i]

    def compute_COM(self, atoms):
        
        xyz = np.zeros(3)
        M = 0

        for i in atoms:

            symbol = self.traj.topology.atom(i).element.symbol
            m = mdlv.element(symbol).atomic_weight
            # this gives positions in nm
            position =  np.array(self.traj.xyz[0][i])

            M += m
            xyz += m * position
            
        xyz *= 1/M
        
        return xyz

    def compute_lamda_i(self, r_i):

        x_i = float((r_i - self.Rmin) / (self.Rmax - self.Rmin))

        lamda_i = -6*((x_i)**5) + 15*((x_i)**4) - 10*((x_i)**3) + 1

        d_lamda_i = -30*(x_i)**4  + 60*(x_i)**3 - 30*(x_i)**2

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


